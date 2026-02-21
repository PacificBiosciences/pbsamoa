#ifndef PBCOPPER_PARALLEL_THREADPOOL_H
#define PBCOPPER_PARALLEL_THREADPOOL_H

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <stop_token>
#include <thread>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Parallel {

///
/// \brief Strong, private type representing a thread index within the thread pool.
///
struct ThreadIndex
{
    ///
    /// \brief Explicit conversion to underlying integer value (0-based thread index).
    ///
    /// \return The thread index as a std::size_t
    ///
    [[nodiscard]] explicit constexpr operator std::size_t() const noexcept { return value_; }

    ///
    /// \brief Get the underlying integer value (0-based thread index).
    ///
    /// \return The thread index as a std::size_t
    ///
    [[nodiscard]] constexpr std::size_t Value() const noexcept { return value_; }

private:
    std::size_t value_;

    ///
    /// \brief Construct a ThreadIndex from an integer value
    ///
    /// \param value The thread index value (must be >= 0)
    ///
    explicit constexpr ThreadIndex(std::size_t value) noexcept : value_{value} {}

    /// Only ThreadPool can create ThreadIndex instances
    template <typename>
    friend class ThreadPool;
};

///
/// \brief Operational metrics for monitoring ThreadPool performance and health
///
/// Provides real-time insights into task throughput, queue utilization, and failure rates. All
/// metrics are updated atomically for thread-safe access without locks.
///
/// **Usage:** Enable metrics collection by enabling Config::EnableMetrics, then periodically call
/// GetSnapshot() to retrieve a consistent snapshot.
///
/// \note Can only be instantiated by ThreadPool.
///
struct Metrics
{
    ///
    /// \brief Immutable, atomic snapshot of ThreadPool metrics at a point in time
    ///
    struct Snapshot
    {
        /// Total tasks submitted via Submit() or SubmitWithPriority()
        std::size_t TasksSubmitted{0};

        /// Total tasks that completed successfully (no exception)
        /// For producer-consumer pools (ResultType != void), this increments when results are
        /// consumed (future::get) and may lag execution completion.
        std::size_t TasksCompleted{0};

        /// Total tasks that threw exceptions during execution
        /// For producer-consumer pools (ResultType != void), this increments when results are
        /// consumed (future::get) and may lag execution completion.
        std::size_t TasksFailed{0};

        /// Current number of tasks in the work queue (head_)
        std::size_t CurrentQueueDepth{0};

        /// Maximum work queue depth observed since pool creation
        std::size_t PeakQueueDepth{0};

        /// Current number of tasks actively executing in worker threads
        std::size_t CurrentActiveTasks{0};

        /// Maximum active tasks observed simultaneously
        std::size_t PeakActiveTasks{0};

        /// Current number of results awaiting consumption (producer-consumer mode only)
        /// Updated when futures are enqueued and consumed; may be non-zero after execution if the
        /// consumer lags.
        /// Always 0 for fire-and-forget mode (ResultType = void)
        std::size_t CurrentResultQueueDepth{0};

        /// Maximum result queue depth observed since pool creation (producer-consumer mode only)
        std::size_t PeakResultQueueDepth{0};
    };

    /// Atomically capture current metrics as an immutable snapshot
    [[nodiscard]] Snapshot GetSnapshot() const noexcept
    {
        return Snapshot{
            .TasksSubmitted = tasksSubmitted.load(std::memory_order_relaxed),
            .TasksCompleted = tasksCompleted.load(std::memory_order_relaxed),
            .TasksFailed = tasksFailed.load(std::memory_order_relaxed),
            .CurrentQueueDepth = currentQueueDepth.load(std::memory_order_relaxed),
            .PeakQueueDepth = peakQueueDepth.load(std::memory_order_relaxed),
            .CurrentActiveTasks = currentActiveTasks.load(std::memory_order_relaxed),
            .PeakActiveTasks = peakActiveTasks.load(std::memory_order_relaxed),
            .CurrentResultQueueDepth = currentResultQueueDepth.load(std::memory_order_relaxed),
            .PeakResultQueueDepth = peakResultQueueDepth.load(std::memory_order_relaxed),
        };
    }

private:
    std::atomic<std::size_t> tasksSubmitted{0};
    std::atomic<std::size_t> tasksCompleted{0};
    std::atomic<std::size_t> tasksFailed{0};
    std::atomic<std::size_t> currentQueueDepth{0};
    std::atomic<std::size_t> peakQueueDepth{0};
    std::atomic<std::size_t> currentActiveTasks{0};
    std::atomic<std::size_t> peakActiveTasks{0};
    std::atomic<std::size_t> currentResultQueueDepth{0};
    std::atomic<std::size_t> peakResultQueueDepth{0};

    template <typename>
    friend class ThreadPool;
};

///
/// \brief High-performance thread pool with support for fire-and-forget and producer-consumer
/// patterns
///
/// ThreadPool is a zero-overhead abstraction for parallel task execution that unifies two
/// fundamental parallelism patterns through a single, type-safe interface controlled by the
/// ResultType template parameter.
///
/// # Operating Modes
///
/// ## Fire-and-Forget Mode (ResultType = void)
///
/// Submit independent tasks that execute asynchronously without result collection.
///
/// **Characteristics:**
/// - Tasks execute as soon as worker threads are available
/// - Completion detected via Finalize() which blocks until all tasks finish
///
/// ## Producer-Consumer Mode (ResultType != void)
///
/// Submit tasks that produce results, which are consumed by a dedicated consumer thread via
/// ConsumeWith().
///
/// **Characteristics:**
/// - Results are dequeued according to scheduling priority (FIFO within each priority level)
/// - Consumer receives results in batches to reduce lock contention
/// - Back-pressure: producers block when result queue fills
/// - Submit() blocks until a consumer thread has entered ConsumeWith()
///
/// # Priority Scheduling
///
/// Tasks submitted via Submit() execute in FIFO order. SubmitWithPriority() enables explicit
/// scheduling control when some work must run sooner than others.
///
/// **Characteristics:**
/// - Lower numeric priority values run before higher values
/// - Equal priorities preserve FIFO ordering (stable priority queue)
/// - Submit() is equivalent to SubmitWithPriority(0, ...)
///
/// # Thread Index Support
///
/// Every task receives a zero-based thread index (0 to NumThreads-1) that identifies which worker
/// thread is executing it. This enables thread-local optimizations without external
/// synchronization. Thread indices are optionally passed to the task as the first parameter.
///
/// # Exception Handling Model
///
/// ThreadPool implements a fail-fast exception policy:
///
/// 1. **First Exception Wins**: The first exception from any source (worker thread, submitted task,
///                              or thread exit callback) is captured
/// 2. **Immediate Shutdown**:   All worker threads are signaled to stop
/// 3. **Propagation**:          Exception is rethrown to the caller of:
///                               - Finalize() for fire-and-forget mode
///                               - Dispatch() for fire-and-forget mode
///                               - ConsumeWith() for producer-consumer mode
///                               - Submit() if pool already failed
/// 4. **Cleanup**:              All in-flight tasks complete, queued tasks may be dropped,
///                              and no new tasks are accepted
///
/// **Important:** Multiple exceptions may occur, but only the first is preserved. Subsequent
///                exceptions are silently discarded. This is by design to avoid exception ordering
///                non-determinism.
///
/// # Synchronization and Thread Safety
///
/// - **Submit()**:      Thread-safe. Multiple threads may submit concurrently.
/// - **ConsumeWith()**: Not thread-safe. Must be called from a single consumer thread.
/// - **Finalize()**:    Thread-safe and idempotent, but must not be called from a worker thread.
/// - **Destructor**:    Automatically calls Finalize() if not already called; must not run on a
///                      worker thread.
///
/// **Producer back-pressure:** Submit() blocks when work queue reaches capacity. In
///                             producer-consumer mode, Submit() also blocks until a consumer
///                             thread has entered ConsumeWith().
///
/// **Consumer back-pressure:** In producer-consumer mode, workers block when result queue reaches
///                             capacity, preventing producers from overwhelming the consumer.
///
/// \tparam ResultType Return type of submitted tasks (void for fire-and-forget mode).
///
template <typename ResultType = void>
class ThreadPool
{
    /// Compile-time mode detection for clearer conditionals
    static constexpr bool IsFireAndForget = std::is_void_v<ResultType>;
    static constexpr bool IsProducerConsumer = !IsFireAndForget;

public:
    ///
    /// \brief Configuration for ThreadPool construction
    ///
    struct Config
    {
        /// \brief Number of worker threads to spawn (default = 1)
        ///
        /// **Must be > 0.** Constructor throws std::invalid_argument if zero.
        ///
        std::size_t NumThreads{1};

        /// \brief Capacity multiplier for internal work and result queues (default = 3).
        ///
        /// **Must be > 0.** Constructor throws std::invalid_argument if zero.
        ///
        std::size_t QueueMultiplier{3};

        /// \brief Optional callback invoked when each worker thread exits
        ///
        /// \note **When called:**
        ///       - After the thread's work loop completes normally
        ///       - Before the thread is joined in Finalize()
        ///       - With the thread's ThreadIndex (0 to NumThreads-1)
        ///
        /// \note **When NOT called (best-effort):**
        ///       - Generally skipped if an exception occurred in any thread (triggers fail-fast shutdown)
        ///       - Generally skipped if the destructor forces fail-fast shutdown
        ///       - **RACE CONDITION:** Due to TOCTOU, the callback may occasionally execute even after
        ///         an exception was captured. The implementation checks exception state before calling,
        ///         but another thread may set an exception between the check and callback invocation.
        ///         If your callback requires strict "only on success" guarantees, implement your own
        ///         exception state tracking.
        ///
        /// \note **Exception behavior:**
        ///       - Exceptions thrown by onThreadExit are captured and treated as task exceptions:
        ///         - Triggers fail-fast shutdown of all workers
        ///         - Exception is rethrown from Finalize()
        ///         - Only the first exception is preserved
        ///
        /// \note **Thread safety:** The callback may be invoked concurrently by multiple threads as
        ///       they exit. Synchronize access to shared state if needed.
        ///
        std::function<void(ThreadIndex threadIndex)> OnThreadExit = {};

        /// \brief Enable operational metrics collection (default = false)
        ///
        /// When enabled, the pool tracks task submission, completion, failure rates, and queue
        /// depths. Metrics can be retrieved via ThreadPool::GetMetrics().
        ///
        bool EnableMetrics{false};
    };

    ///
    /// \brief Construct and start a ThreadPool with the specified configuration
    ///
    /// \param config Configuration specifying thread count, queue capacity, and optional thread
    ///               exit callback. See ThreadPool::Config for detailed parameter documentation.
    ///
    /// \throws std::invalid_argument if config.NumThreads == 0 or config.QueueMultiplier == 0
    /// \throws std::system_error if thread creation fails
    ///
    /// \note **Construction behavior:**
    ///       - Validates configuration parameters
    ///       - Spawns numThreads worker threads immediately
    ///       - Worker threads enter their event loop, waiting for tasks
    ///       - Pool is ready to accept Submit() calls immediately after construction
    ///
    explicit ThreadPool(Config config = {});

    ///
    /// \brief Destructor - ensures all work completes and threads are joined
    ///
    /// \warning If an exception occurred in a worker thread and you did not call Finalize()
    ///          explicitly, the exception will be silently suppressed by the destructor.
    ///          This may hide bugs.
    /// \warning Destroying a ThreadPool from one of its worker threads is undefined behavior and
    ///          can lead to use-after-free. Always destroy from a non-worker thread.
    ///
    /// \note **Automatic finalization:**
    ///       If Finalize() has not been explicitly called, destructor automatically invokes it to
    ///       ensure:
    ///       - All in-flight tasks complete execution (queued tasks may be dropped on fail-fast)
    ///       - All worker threads are properly joined
    ///
    /// \note **Thread safety:**
    ///       The destructor is not thread-safe. Destroying a ThreadPool while other threads might
    ///       be calling Submit() or ConsumeWith() is undefined behavior. It must also not be
    ///       called from a worker thread.
    ///
    /// \note **Performance note:**
    ///       Destruction blocks until all tasks complete. For pools with long-running tasks,
    ///       destruction may take significant time. Ensure the destroying thread can afford to wait.
    ///
    ~ThreadPool();

    /// Non-copyable, non-movable
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool& operator=(ThreadPool&&) = delete;

    ///
    /// \brief Submit a task for execution by the thread pool
    ///
    /// Automatically detects whether the callable accepts a ThreadIndex parameter. If it does, the
    /// ThreadIndex will be passed as the first argument. Otherwise, the callable is wrapped to
    /// accept but ignore the ThreadIndex parameter.
    ///
    /// \tparam F Callable type (lambda, function pointer, functor, std::function)
    /// \tparam Args Argument types to forward to the callable
    ///
    /// \param func Callable with one of these signatures:
    ///             - ResultType(ThreadIndex threadIndex, Args...) - uses ThreadIndex
    ///             - ResultType(Args...) - ignores ThreadIndex
    ///               Must be movable or copyable. Captures are moved/copied once.
    ///               For large captures, consider capturing by pointer/reference.
    /// \note **ThreadIndex detection:**
    ///       The callable must explicitly accept ThreadIndex to receive it.
    /// \param args Arguments perfect-forwarded to func (after the optional ThreadIndex).
    ///             Each argument is moved/copied once into the internal task storage.
    ///
    /// \throws std::runtime_error if Finalize() was already called
    /// \throws The first exception captured by the pool if a prior task or callback failed
    ///         (fail-fast: pool stops accepting work after first error)
    /// \throws std::runtime_error if ResultType != void and the pool is shutting down before a
    ///                            consumer thread has entered ConsumeWith()
    ///
    /// \note **Execution model:**
    ///       - **Fire-and-forget mode:**   Task executes asynchronously.
    ///                                     Completion is detected via Finalize().
    ///       - **Producer-consumer mode:** Task's return value is queued for consumption via
    ///                                     ConsumeWith() in scheduler order (priority, FIFO
    ///                                     within priority).
    ///
    /// \note **Blocking behavior:**
    ///       Submit() blocks if the work queue is full (NumThreads * QueueMultiplier tasks already
    ///       queued). This provides automatic back-pressure. The call unblocks as soon as a worker
    ///       dequeues a task.
    ///
    /// \note **Producer-consumer guard:**
    ///       For ResultType != void, Submit() blocks until a consumer thread has entered
    ///       ConsumeWith(), to prevent deadlock when the result queue fills. If the pool is
    ///       shutting down, Submit() throws.
    ///
    /// \note **Exception capture:**
    ///       Exceptions thrown by func are captured and stored. The first exception stops the pool
    ///       and is rethrown from Finalize() or ConsumeWith().
    ///       Submit() may also rethrow if called after failure.
    ///
    /// \note **Thread safety:**
    ///       Multiple threads may call Submit() concurrently.
    ///       Internal synchronization ensures correct queuing.
    ///
    /// \note **Priority:**
    ///       Equivalent to SubmitWithPriority(0, ...).
    ///       Use SubmitWithPriority() to provide explicit priorities.
    ///
    template <typename F, typename... Args>
    void Submit(F&& func, Args&&... args);

    ///
    /// \brief Submit a task with an explicit scheduling priority
    ///
    /// \tparam F Callable type (lambda, function pointer, functor, std::function)
    /// \tparam Args Argument types to forward to the callable
    ///
    /// \param priority Task priority. Lower values run sooner. Equal priorities
    ///                 preserve FIFO ordering. Priority does not affect failure semantics or queue
    ///                 capacity.
    /// \param func Callable as described in Submit() documentation.
    /// \param args Arguments perfect-forwarded to func (after the optional ThreadIndex).
    ///
    /// \note **Priority:**
    ///       Lower priority values execute before higher values. Tasks submitted with identical
    ///       priorities run in submission order.
    ///
    /// \note **Behavior:**
    ///       Same blocking and exception semantics as Submit(), including the producer-consumer
    ///       guard.
    ///
    template <typename F, typename... Args>
    void SubmitWithPriority(std::int32_t priority, F&& func, Args&&... args);

    ///
    /// \brief Dispatch a fixed number of tasks with exception handling and synchronization
    ///
    /// \param numEntries Number of tasks to dispatch (must be >= 0)
    /// \param callback Function to call for each task index: void(std::int32_t index)
    ///
    /// \throws std::invalid_argument if numEntries < 0
    /// \throws Any exception thrown by the callback (first exception wins)
    /// \throws std::runtime_error if Finalize() was already called
    /// \throws The first exception captured by the pool if a prior task or callback failed
    /// \throws std::runtime_error if called on producer-consumer pool (ResultType != void)
    ///
    /// \note **Dispatch functionality:** This method provides dispatch functionality for
    ///       fire-and-forget style task execution. It submits tasks for each index
    ///       from 0 to numEntries-1, waits for all to complete, and handles exceptions with a
    ///       "first exception wins" policy. If an exception occurs, Dispatch waits for in-flight
    ///       tasks to finish, and queued tasks may be abandoned due to fail-fast shutdown.
    ///
    /// \note **Priority:**
    ///       Equivalent to DispatchWithPriority(numEntries, 0, ...).
    ///
    /// \note **Availability:**
    ///       This is only available for fire-and-forget mode (ResultType = void).
    ///
    template <typename F>
        requires std::is_void_v<ResultType>
    void Dispatch(std::int32_t numEntries, F&& callback);

    ///
    /// \brief Dispatch tasks with an explicit scheduling priority
    ///
    /// Identical to Dispatch() but tasks are submitted via SubmitWithPriority() using the provided
    /// priority value.
    ///
    /// \param numEntries Number of tasks to dispatch (must be >= 0)
    /// \param priority Priority passed to SubmitWithPriority(). Lower values run sooner.
    /// \param callback Function to call for each task index: void(std::int32_t index)
    ///
    /// \throws std::invalid_argument if numEntries < 0
    /// \throws Any exception thrown by the callback (first exception wins)
    /// \throws std::runtime_error if Finalize() was already called
    /// \throws The first exception captured by the pool if a prior task or callback failed
    /// \throws std::runtime_error if called on producer-consumer pool (ResultType != void)
    ///
    /// \note **Availability:** This is only available for fire-and-forget mode (ResultType = void).
    ///
    template <typename F>
        requires std::is_void_v<ResultType>
    void DispatchWithPriority(std::int32_t numEntries, std::int32_t priority, F&& callback);

    ///
    /// \brief Finalize the pool: complete all work, join threads, propagate exceptions
    ///
    /// \throws Any exception type that was thrown by:
    ///         - A submitted task (via Submit())
    ///         - A thread exit callback (Config::onThreadExit)
    ///         - Internal worker thread logic
    /// \throws std::runtime_error if called from a worker thread
    ///
    /// \note **Idempotency:**
    ///       Finalize() is safe to call multiple times. The first call performs all work;
    ///       subsequent calls return immediately but still rethrow any captured exception.
    ///
    /// \note **Blocking behavior:**
    ///       This is a synchronous blocking call that may take arbitrarily long depending on task
    ///       execution time. Ensure the calling thread can afford to wait.
    ///
    /// \note **Mode-specific behavior:**
    ///       - Fire-and-forget mode (ResultType=void):
    ///         - Blocks until all tasks complete
    ///         - If an exception occurs, pending queued tasks are cleared (fail-fast)
    ///         - No consumer coordination needed
    ///       - Producer-consumer mode (ResultType!=void):
    ///         - Blocks until all tasks complete and results are queued
    ///         - If an exception occurs, pending queued tasks are cleared (fail-fast)
    ///           and results for those tasks are not produced
    ///         - Sends shutdown signal (empty optional) to result queue
    ///         - Consumer's next ConsumeWith() call will return false
    ///         - Consumer thread should check return value and exit
    ///
    /// \note **Exception handling:**
    ///       If any exception occurred (in any worker thread, submitted task, or thread exit
    ///       callback), Finalize() rethrows the **first** exception that was captured.
    ///       Subsequent exceptions are silently discarded to avoid non-determinism.
    ///
    /// \note **Best practice:**
    ///       Always call Finalize() explicitly (not just in destructor) to handle exceptions properly.
    ///       The destructor suppresses exceptions.
    ///
    /// \note **Thread safety:**
    ///       Finalize() is thread-safe and can be called from any non-worker thread.
    ///       However, calling Finalize() concurrently with Submit() may cause Submit() to throw if
    ///       the pool shuts down between the calls.
    ///
    /// \note **Post-finalization:**
    ///       After Finalize() completes, Submit() will throw std::runtime_error.
    ///       The pool cannot be restarted—create a new instance.
    ///
    void Finalize();

    ///
    /// \brief Consume a batch of results from the result queue (producer-consumer mode only)
    ///
    /// \warning **CRITICAL: Must Call Before Submitting Tasks**
    ///          In producer-consumer mode, you MUST start a dedicated thread calling ConsumeWith()
    ///          in a loop BEFORE submitting any tasks. If the consumer is not actively running,
    ///          the pool will deadlock once the result queue fills (NumThreads * QueueMultiplier).
    ///          Submit() enforces this by blocking until a consumer thread has entered
    ///          ConsumeWith(). If you need strict ordering, use a startup handshake to ensure the
    ///          consumer thread has started. See class-level documentation for examples.
    ///
    /// \tparam Consumer Callable type that processes each result
    ///
    /// \param consumer Callable with signature: void(ResultType) or void(const
    ///                 ResultType&). Supports both lvalue references (named
    ///                 callables like lambdas with captures, functors, or
    ///                 function objects) and rvalue references (temporary
    ///                 callables like std::function). Automatic lifetime
    ///                 management ensures proper handling of both reference types.
    ///                 Invoked once per result in scheduler order. Should be
    ///                 efficient to maintain throughput. For expensive
    ///                 processing, consider queueing results to another thread
    ///                 pool.
    ///
    /// \returns true if results were consumed and more may arrive (keep calling),
    ///          false if Finalize() was called and result queue is drained (stop
    ///          consuming)
    ///
    /// \throws Any exception thrown by a submitted task (via std::future::get())
    /// \throws Any exception thrown by the consumer callable
    /// \throws The first exception captured by the pool if an exception occurred elsewhere
    ///         (via ThrowIfAborted())
    ///
    /// \note **Only available when ResultType != void.** Attempting to call this on a
    ///       fire-and-forget pool (ResultType=void) will result in a compile-time error
    ///       due to the requires clause.
    ///
    /// \note **Consumption model:**
    ///       - **Blocks** until at least one result is available
    ///       - **Atomically swaps** the entire result queue into local storage
    ///       - **Invokes consumer** once for each result in FIFO order
    ///       - **Returns** true if more results may arrive, false if pool finalized
    ///
    /// \note **Batching behavior:**
    ///       ConsumeWith() processes ALL available results in a single call, not just one. This
    ///       batching reduces lock contention by minimizing mutex acquisitions.
    ///       For a high-throughput pipeline, this can significantly improve performance.
    ///
    /// \note **Typical usage pattern:**
    ///       Create a dedicated consumer thread that calls ConsumeWith() in a loop until it returns
    ///       false, indicating the pool has been finalized.
    ///
    /// \note **Execution order:**
    ///       Results are consumed in the order tasks leave the scheduler: lower priority values
    ///       first, FIFO within the same priority. Batching means multiple results are processed
    ///       before the next ConsumeWith() call.
    ///
    /// \note **Performance considerations:**
    ///       - Consumer callable is invoked OUTSIDE the internal lock
    ///       - Each result must be unwrapped from a std::future (may block briefly)
    ///       - Expensive consumer processing does NOT block producers
    ///       - For very high throughput, ensure consumer can keep pace with producers
    ///
    /// \note **Exception handling:**
    ///       - If a task threw an exception, future::get() rethrows it from ConsumeWith()
    ///       - If the consumer throws, the exception propagates to the caller
    ///       - In both cases, the pool enters fail-fast mode and no new work is accepted
    ///
    /// \note **Thread safety:**
    ///       ConsumeWith() is NOT thread-safe. It must be called from a single dedicated consumer
    ///       thread. Calling from multiple threads concurrently will cause data races and undefined
    ///       behavior.
    ///
    /// \note **Compiler enforcement:**
    ///       This method is disabled at compile-time for fire-and-forget pools (ResultType=void) via
    ///       the requires clause. Attempting to call it will produce a clear compiler error.
    ///
    template <typename Consumer>
        requires std::invocable<Consumer, ResultType> && (!std::is_void_v<ResultType>)
    bool ConsumeWith(Consumer&& consumer);

    ///
    /// \brief Retrieve a snapshot of operational metrics
    /// \returns Metrics::Snapshot containing current metric values
    ///
    /// \note **Thread safety:** Safe to call from any thread at any time
    ///
    /// \note **Performance:** Lock-free atomic loads, ~50-100ns total latency
    ///
    /// \note **Availability:**
    ///       Returns valid data only if Config::EnableMetrics was set to true during construction.
    ///       Otherwise, all values in the snapshot will be zero.
    ///
    /// \note **Producer-consumer semantics:**
    ///       For ResultType != void, TasksCompleted/TasksFailed are recorded when results are
    ///       consumed (future::get), not at task execution completion.
    ///
    [[nodiscard]] Metrics::Snapshot GetMetrics() const noexcept
    {
        if (metrics_) {
            return metrics_->GetSnapshot();
        }
        return Metrics::Snapshot{};
    }

private:
    using TaskType = std::packaged_task<ResultType(ThreadIndex)>;
    using OptionalTask = std::optional<TaskType>;
    using OptionalFuture =
        std::conditional_t<IsProducerConsumer, std::optional<std::future<ResultType>>,
                           std::monostate>;

    /// A prioritized task is a task with a priority and a sequence number.
    /// The priority is used to determine the order in which the task is executed.
    /// The sequence number is used to determine the order in which the task is executed
    /// if the priorities are the same.
    struct PrioritizedTask
    {
        std::int32_t Priority;
        std::size_t Sequence;
        TaskType Task;
    };

    struct TaskPriorityCompare
    {
        [[nodiscard]] bool operator()(const PrioritizedTask& lhs,
                                      const PrioritizedTask& rhs) const noexcept
        {
            if (lhs.Priority != rhs.Priority) {
                return lhs.Priority > rhs.Priority;
            }
            return lhs.Sequence > rhs.Sequence;
        }
    };

    template <typename F, typename... Args>
    [[nodiscard]] static TaskType WrapCallable(F&& func, Args&&... args)
    {
        static_assert(
            std::is_invocable_v<F&, ThreadIndex, Args...> || std::is_invocable_v<F&, Args...>,
            "ThreadPool::Submit requires a callable invocable with (ThreadIndex, Args...) "
            "or (Args...).");

        if constexpr (std::is_invocable_v<F&, ThreadIndex, Args...>) {
            return TaskType{
                [f = std::forward<F>(func), ... capturedArgs = std::forward<Args>(args)](
                    ThreadIndex threadIndex) mutable -> decltype(auto) {
                    return f(threadIndex, std::move(capturedArgs)...);
                }};
        } else {
            return TaskType{
                [f = std::forward<F>(func), ... capturedArgs = std::forward<Args>(args)](
                    [[maybe_unused]] ThreadIndex threadIndex) mutable -> decltype(auto) {
                    return f(std::move(capturedArgs)...);
                }};
        }
    }

    /// Worker threads executing tasks from the work queue.
    std::vector<std::jthread> threads_;
    /// Priority queue of pending tasks (min-heap by priority, FIFO within same priority).
    std::vector<PrioritizedTask> head_;

    /// Result queue for producer-consumer mode (ResultType != void).
    /// Zero-size (std::monostate) in fire-and-forget mode (ResultType = void).
    /// [[no_unique_address]] eliminates storage overhead when unused.
    [[no_unique_address]] std::conditional_t<IsProducerConsumer, std::deque<OptionalFuture>,
                                             std::monostate> tail_;

    /// Consumer thread ID for single-consumer enforcement (producer-consumer mode only).
    /// Zero-size (std::monostate) in fire-and-forget mode (ResultType = void).
    /// [[no_unique_address]] eliminates storage overhead when unused.
    [[no_unique_address]] std::conditional_t<IsProducerConsumer, std::optional<std::thread::id>,
                                             std::monostate> consumerThreadId_;
    /// Tracks whether a consumer is actively inside ConsumeWith (producer-consumer mode only).
    /// Zero-size (std::monostate) in fire-and-forget mode (ResultType = void).
    [[no_unique_address]] std::conditional_t<IsProducerConsumer, bool, std::monostate>
        consumerActive_{};

    /// Mutex protecting all shared pool state (head_, tail_, threads_, etc.). Mutable to allow
    /// locking in const member functions like IsWorkerThread().
    mutable std::mutex poolMutex_;
    /// Condition variable for notifying worker threads that new work is available or shutdown
    /// is requested.
    std::condition_variable pushed_;
    /// Condition variable for notifying producers (Submit callers) that the input queue has
    /// available capacity.
    std::condition_variable cvInputQueueNotFull_;
    /// Condition variable for notifying the consumer thread (ConsumeWith caller) that results
    /// are available or shutdown is requested.
    std::condition_variable cvResultQueueNotEmpty_;
    /// Condition variable for notifying producers (Submit callers) that a consumer thread has
    /// entered ConsumeWith().
    std::condition_variable cvConsumerRegistered_;
    /// Condition variable for notifying Finalize() that all tasks have completed or the work
    /// queue is drained.
    std::condition_variable cvAllTasksComplete_;
    /// Stop source for signaling worker threads to stop (fail-fast shutdown on exception).
    std::stop_source stopSource_;
    /// Atomic flag indicating if the pool is accepting new job submissions.
    std::atomic<bool> acceptingJobs_{true};
    /// Atomic sequence number for stable priority queue ordering (monotonically increasing).
    std::atomic<std::size_t> nextSequence_{0};
    /// Atomic count of tasks currently executing in worker threads.
    std::atomic<std::size_t> activeTasks_{0};
    /// Maximum capacity of the input work queue (NumThreads * QueueMultiplier).
    std::size_t maxQueueDepth_{0};
    /// Optional callback invoked when each worker thread exits (see Config::OnThreadExit).
    std::function<void(ThreadIndex)> onThreadExit_;
    /// Optional metrics tracker for monitoring pool performance (nullptr if disabled).
    std::unique_ptr<Metrics> metrics_;
    /// First captured exception (fail-fast: only the first exception is preserved).
    std::exception_ptr exc_;
    /// Mutex protecting exc_ and exceptionOccurred_.
    mutable std::mutex exceptionMutex_;
    /// Fast check flag indicating if an exception has occurred (avoids dereferencing exc_).
    bool exceptionOccurred_{false};
    /// Ensures Finalize() logic runs exactly once, even with concurrent calls.
    std::once_flag finalizeOnce_;

    void WorkerLoop(ThreadIndex threadIndex, std::stop_token stopToken);
    OptionalTask PopTask(std::stop_token stopToken);
    void ExecuteTask(TaskType& task, ThreadIndex threadIndex);
    void ProduceTask(TaskType&& task, std::int32_t priority);
    template <typename SubmitFunc, typename F>
        requires std::is_void_v<ResultType>
    void DispatchImpl(std::int32_t numEntries, SubmitFunc submitFunc, F&& callback);
    void SetFirstException(std::exception_ptr ptr = std::current_exception()) noexcept;
    void ThrowIfAborted();
    void SignalNoMoreTasks() noexcept;
    void JoinAllThreads(bool requestStop = false) noexcept;
    void WakeAllWaiters() noexcept;
    template <typename Predicate>
    void WaitForCondition(Predicate&& pred);

    /// \brief Helper: update metrics atomically (no-op if metrics disabled)
    void UpdateMetricsOnSubmit() noexcept
    {
        if (metrics_) {
            metrics_->tasksSubmitted.fetch_add(1, std::memory_order_relaxed);
        }
    }

    /// \brief Helper: update a peak metric using CAS loop (no-op if metrics disabled)
    void UpdatePeakMetric(std::atomic<std::size_t>& peak, std::size_t newValue) noexcept
    {
        std::size_t current{peak.load(std::memory_order_relaxed)};
        while (newValue > current &&
               !peak.compare_exchange_weak(current, newValue, std::memory_order_relaxed)) {
        }
    }

    void UpdateMetricsOnQueueChange(std::size_t newDepth) noexcept
    {
        if (metrics_) {
            metrics_->currentQueueDepth.store(newDepth, std::memory_order_relaxed);
            UpdatePeakMetric(metrics_->peakQueueDepth, newDepth);
        }
    }

    void IncrementResultQueueDepth() noexcept
    {
        if (metrics_) {
            const std::size_t newDepth{
                metrics_->currentResultQueueDepth.fetch_add(1, std::memory_order_relaxed) + 1};
            UpdatePeakMetric(metrics_->peakResultQueueDepth, newDepth);
        }
    }

    void DecrementResultQueueDepth() noexcept
    {
        if (metrics_) {
            metrics_->currentResultQueueDepth.fetch_sub(1, std::memory_order_relaxed);
        }
    }

    void UpdateMetricsOnTaskStart() noexcept
    {
        if (metrics_) {
            const std::size_t active{activeTasks_.load(std::memory_order_relaxed)};
            metrics_->currentActiveTasks.store(active, std::memory_order_relaxed);
            UpdatePeakMetric(metrics_->peakActiveTasks, active);
        }
    }

    void UpdateMetricsOnTaskComplete(std::optional<bool> failed = std::nullopt) noexcept
    {
        if (metrics_) {
            if (failed.has_value()) {
                (*failed ? metrics_->tasksFailed : metrics_->tasksCompleted)
                    .fetch_add(1, std::memory_order_relaxed);
            }
            metrics_->currentActiveTasks.store(activeTasks_.load(std::memory_order_relaxed),
                                               std::memory_order_relaxed);
        }
    }

    /// \brief Helper: detect if the caller is one of our worker threads
    ///
    /// \return true if the caller is one of our worker threads, false otherwise
    ///
    [[nodiscard]] bool IsWorkerThread() const
    {
        std::lock_guard lock{poolMutex_};
        const std::thread::id self{std::this_thread::get_id()};
        for (const auto& thread : threads_) {
            if (thread.get_id() == self) {
                return true;
            }
        }
        return false;
    }
};

// ========================================================================
// Implementation
// ========================================================================

template <typename ResultType>
ThreadPool<ResultType>::ThreadPool(Config config)
    : onThreadExit_{std::move(config.OnThreadExit)}
    , metrics_{config.EnableMetrics ? std::make_unique<Metrics>() : nullptr}
{
    if (config.NumThreads == 0) {
        throw std::invalid_argument{"ThreadPool: numThreads must be > 0"};
    }
    if (config.QueueMultiplier == 0) {
        throw std::invalid_argument{"ThreadPool: queueMultiplier must be > 0"};
    }
    if (config.NumThreads > (std::numeric_limits<std::size_t>::max() / config.QueueMultiplier)) {
        throw std::invalid_argument{
            "ThreadPool: numThreads * queueMultiplier would overflow size_t"};
    }

    maxQueueDepth_ = config.NumThreads * config.QueueMultiplier;
    threads_.reserve(config.NumThreads);
    head_.reserve(maxQueueDepth_);

    try {
        for (std::size_t i{0}; i < config.NumThreads; ++i) {
            threads_.emplace_back([this, threadIndex{i}](std::stop_token token) {
                WorkerLoop(ThreadIndex{threadIndex}, token);
            });
        }
    } catch (...) {
        // Ensure running threads are woken up so they can exit and be joined
        stopSource_.request_stop();
        {
            std::lock_guard lock{poolMutex_};
            pushed_.notify_all();  // Wake workers so they see the stop token
        }
        throw;  // Vector destructor handles joining
    }
}

template <typename ResultType>
ThreadPool<ResultType>::~ThreadPool()
{
    try {
        if (acceptingJobs_.load(std::memory_order_acquire)) {
            Finalize();
        }
    } catch (...) {
    }
}

template <typename ResultType>
template <typename F, typename... Args>
void ThreadPool<ResultType>::Submit(F&& func, Args&&... args)
{
    SubmitWithPriority(0, std::forward<F>(func), std::forward<Args>(args)...);
}

template <typename ResultType>
template <typename F, typename... Args>
void ThreadPool<ResultType>::SubmitWithPriority(std::int32_t priority, F&& func, Args&&... args)
{
    TaskType task{WrapCallable(std::forward<F>(func), std::forward<Args>(args)...)};
    ProduceTask(std::move(task), priority);
}

template <typename ResultType>
void ThreadPool<ResultType>::Finalize()
{
    // Disallow invoking Finalize() from a worker thread to avoid self-deadlock.
    if (IsWorkerThread()) {
        throw std::runtime_error{"ThreadPool::Finalize() may not be called from a worker thread"};
    }

    // Ensure finalization logic runs exactly once, even with concurrent calls
    std::call_once(finalizeOnce_, [this]() {
        SignalNoMoreTasks();

        // Fire-and-forget mode: drain work queue before waiting for active tasks
        if constexpr (IsFireAndForget) {
            WaitForCondition(
                [this] { return std::ranges::empty(head_) || stopSource_.stop_requested(); });

            if (stopSource_.stop_requested()) {
                std::lock_guard lock{poolMutex_};
                head_.clear();
                cvInputQueueNotFull_.notify_all();
            }
        }

        // Wait for all in-flight tasks to complete
        WaitForCondition([this] { return activeTasks_.load(std::memory_order_acquire) == 0; });

        JoinAllThreads();

        // Ensure metrics reflect a fully idle pool after finalization
        UpdateMetricsOnTaskComplete();

        // Producer-consumer mode: signal consumer that no more results will arrive
        if constexpr (IsProducerConsumer) {
            {
                std::lock_guard lock{poolMutex_};
                tail_.push_back(std::nullopt);
                cvResultQueueNotEmpty_.notify_all();
            }
        }
    });

    // Always check for exceptions, even on subsequent calls
    ThrowIfAborted();
}

template <typename ResultType>
template <typename F>
    requires std::is_void_v<ResultType>
void ThreadPool<ResultType>::Dispatch(std::int32_t numEntries, F&& callback)
{
    const auto submitTask = [this](auto&& task) {
        this->Submit(std::forward<decltype(task)>(task));
    };
    DispatchImpl(numEntries, submitTask, std::forward<F>(callback));
}

template <typename ResultType>
template <typename F>
    requires std::is_void_v<ResultType>
void ThreadPool<ResultType>::DispatchWithPriority(std::int32_t numEntries, std::int32_t priority,
                                                  F&& callback)
{
    const auto submitTask = [this, priority](auto&& task) {
        this->SubmitWithPriority(priority, std::forward<decltype(task)>(task));
    };
    DispatchImpl(numEntries, submitTask, std::forward<F>(callback));
}

template <typename ResultType>
template <typename SubmitFunc, typename F>
    requires std::is_void_v<ResultType>
void ThreadPool<ResultType>::DispatchImpl(std::int32_t numEntries, SubmitFunc submitFunc,
                                          F&& callback)
{
    if (numEntries < 0) {
        throw std::invalid_argument{"ThreadPool::Dispatch: numEntries must be >= 0"};
    }
    if (numEntries == 0) {
        return;
    }

    auto callbackHolder = std::make_shared<std::decay_t<F>>(std::forward<F>(callback));

    /// Shared state for coordinating a batch of dispatched tasks
    ///
    /// This structure enables Dispatch() to submit multiple independent tasks and wait for all
    /// to complete while implementing "first exception wins" semantics. It is shared across all
    /// tasks in the batch via shared_ptr.
    ///
    /// **Responsibilities:**
    /// - Track completion: Outstanding decrements as each task finishes
    /// - Capture exceptions: ExceptionOnce ensures only the first exception is stored
    /// - Synchronization: Mutex/Cv allow the calling thread to block until Outstanding reaches zero
    ///
    /// **Lifetime:** Kept alive by both the calling thread and all submitted task lambdas
    struct DispatchState
    {
        explicit DispatchState(std::int32_t count) : Outstanding{count} {}

        DispatchState(const DispatchState&) = delete;
        DispatchState& operator=(const DispatchState&) = delete;
        DispatchState(DispatchState&&) = delete;
        DispatchState& operator=(DispatchState&&) = delete;

        /// Number of tasks not yet completed (decremented atomically as tasks finish)
        std::atomic<std::int32_t> Outstanding;
        /// Ensures only the first exception is captured (subsequent exceptions are discarded)
        std::once_flag ExceptionOnce;
        /// The first exception thrown by any task in this dispatch batch
        std::exception_ptr Exception;
        /// Protects access to the condition variable
        std::mutex Mutex;
        /// Notified when Outstanding reaches zero, waking the thread waiting in Dispatch()
        std::condition_variable Cv;
    };

    auto state{std::make_shared<DispatchState>(numEntries)};

    struct DispatchCompletion
    {
        explicit DispatchCompletion(std::shared_ptr<DispatchState> state) noexcept
            : State{std::move(state)}
        {
        }

        DispatchCompletion(const DispatchCompletion&) = delete;
        DispatchCompletion& operator=(const DispatchCompletion&) = delete;

        DispatchCompletion(DispatchCompletion&& other) noexcept
            : State{std::move(other.State)}, Completed{other.Completed}
        {
            other.Completed = true;
        }

        DispatchCompletion& operator=(DispatchCompletion&& other) noexcept
        {
            if (this != &other) {
                Complete();
                State = std::move(other.State);
                Completed = other.Completed;
                other.Completed = true;
            }
            return *this;
        }

        void Complete() noexcept
        {
            if (Completed) {
                return;
            }
            Completed = true;
            if (State && (State->Outstanding.fetch_sub(1, std::memory_order_acq_rel) == 1)) {
                std::lock_guard lock{State->Mutex};
                State->Cv.notify_all();
            }
        }

        ~DispatchCompletion() { Complete(); }

    private:
        std::shared_ptr<DispatchState> State;
        bool Completed{false};
    };

    std::stop_callback stopNotify{stopSource_.get_token(),
                                  [weak = std::weak_ptr<DispatchState>{state}]() noexcept {
                                      if (std::shared_ptr<DispatchState> st = weak.lock()) {
                                          std::lock_guard lock{st->Mutex};
                                          st->Cv.notify_all();
                                      }
                                  }};

    // Helper: capture first exception and trigger fail-fast
    const auto captureException = [this, state]() noexcept {
        std::call_once(state->ExceptionOnce, [state]() {
            state->Exception = std::current_exception();
            if (!state->Exception) {
                state->Exception =
                    std::make_exception_ptr(std::runtime_error{"ThreadPool: dispatch aborted"});
            }
        });
        SetFirstException();  // Triggers stop, which wakes waiter via stop_callback
    };

    // Submit all tasks, capturing exceptions from both submission and execution
    for (std::int32_t i{0}; i < numEntries; ++i) {
        DispatchCompletion completion{state};
        try {
            submitFunc([callbackHolder, captureException, completion = std::move(completion), i]() {
                (void)completion;
                try {
                    std::invoke(*callbackHolder, i);
                } catch (...) {
                    captureException();
                    throw;
                }
            });
        } catch (...) {
            captureException();
        }
    }

    // Wait for all dispatched tasks to complete (or be canceled on fail-fast)
    {
        std::unique_lock lock{state->Mutex};
        state->Cv.wait(lock,
                       [&] { return state->Outstanding.load(std::memory_order_acquire) == 0; });
    }

    // Surface any pool-level exceptions that triggered fail-fast
    ThrowIfAborted();

    if (state->Exception) {
        std::rethrow_exception(state->Exception);
    }
}

template <typename ResultType>
template <typename Consumer>
    requires std::invocable<Consumer, ResultType> && (!std::is_void_v<ResultType>)
bool ThreadPool<ResultType>::ConsumeWith(Consumer&& consumer)
{
    constexpr bool isLvalueConsumer{std::is_lvalue_reference_v<Consumer>};
    using ConsumerStorage =
        std::conditional_t<isLvalueConsumer,
                           std::reference_wrapper<std::remove_reference_t<Consumer>>,
                           std::decay_t<Consumer>>;
    ConsumerStorage consumerCallable{std::forward<Consumer>(consumer)};

    std::deque<OptionalFuture> localTail;

    // Track exception counts to detect stack unwinding reliably.
    const int uncaughtOnEntry{std::uncaught_exceptions()};
    bool isRegisteredConsumer{false};
    bool claimedBatch{false};
    bool consumedShutdownSignal{false};

    // RAII guard: if consumer exits abnormally (throws) or returns with
    // unconsumed results, trigger fail-fast to prevent producer deadlock
    const auto failFastOnExit = [this, &localTail, &isRegisteredConsumer, &claimedBatch,
                                 &consumedShutdownSignal, uncaughtOnEntry](void*) noexcept {
        if (!isRegisteredConsumer) {
            return;  // No consumer was active; nothing to do
        }

        if constexpr (IsProducerConsumer) {
            {
                std::lock_guard lock{poolMutex_};
                consumerActive_ = false;
            }
            cvConsumerRegistered_.notify_all();
        }

        const bool unwinding{std::uncaught_exceptions() > uncaughtOnEntry};
        if (unwinding && claimedBatch) {
            if (std::exception_ptr current{std::current_exception()}) {
                SetFirstException(std::move(current));
            } else {
                SetFirstException(std::make_exception_ptr(
                    std::runtime_error{"ThreadPool: consumer thread unwound"}));
            }
            return;
        }

        if (consumedShutdownSignal) {
            return;  // Normal shutdown
        }

        if (!std::ranges::empty(localTail)) {
            SetFirstException(std::make_exception_ptr(
                std::runtime_error{"ThreadPool: consumer exited with unconsumed results"}));
        }
    };
    std::unique_ptr<void, decltype(failFastOnExit)> guard{this, failFastOnExit};

    {
        std::unique_lock lock{poolMutex_};

        // Enforce single-threaded consumer access
        if constexpr (IsProducerConsumer) {
            const std::thread::id currentThreadId{std::this_thread::get_id()};
            if (!consumerThreadId_) {
                consumerThreadId_ = currentThreadId;
                cvConsumerRegistered_.notify_all();
            } else if (*consumerThreadId_ != currentThreadId) {
                throw std::runtime_error{
                    "ThreadPool::ConsumeWith() must be called from a single thread"};
            }
            isRegisteredConsumer = (*consumerThreadId_ == currentThreadId);
            if (isRegisteredConsumer) {
                consumerActive_ = true;
                cvConsumerRegistered_.notify_all();
            }
        }

        // Predicate must not throw - check conditions and return true to wake up
        cvResultQueueNotEmpty_.wait(lock, [this] {
            // Return true if we should wake up (results available or shutdown)
            return (stopSource_.stop_requested()) || (!std::ranges::empty(tail_));
        });

        // Check for abort after wait returns (predicate must not throw)
        ThrowIfAborted();

        // Batch transfer: swap entire queue to minimize lock contention
        localTail.swap(tail_);
        claimedBatch = true;
        pushed_.notify_all();
    }

    // Process results outside lock
    for (auto& optFuture : localTail) {
        if (!optFuture) {
            localTail.clear();  // Consumed shutdown signal
            claimedBatch = false;
            consumedShutdownSignal = true;
            return false;  // Finalize() sent shutdown signal
        }

        bool decremented{false};
        bool outcomeRecorded{false};
        try {
            decltype(auto) result{optFuture->get()};
            DecrementResultQueueDepth();
            decremented = true;
            UpdateMetricsOnTaskComplete(false);
            outcomeRecorded = true;

            if constexpr (isLvalueConsumer) {
                std::invoke(consumerCallable.get(), std::forward<ResultType>(result));
            } else {
                std::invoke(consumerCallable, std::forward<ResultType>(result));
            }
        } catch (...) {
            if (!decremented) {
                DecrementResultQueueDepth();
            }
            if (!outcomeRecorded) {
                UpdateMetricsOnTaskComplete(true);
            }
            throw;
        }
    }

    claimedBatch = false;
    localTail.clear();  // All consumed successfully
    return true;
}

template <typename ResultType>
void ThreadPool<ResultType>::WorkerLoop(ThreadIndex threadIndex, std::stop_token stopToken)
{
    try {
        while (!stopToken.stop_requested() && !stopSource_.stop_requested()) {
            OptionalTask task{PopTask(stopToken)};
            if (!task) {
                break;
            }
            ExecuteTask(*task, threadIndex);
        }
    } catch (...) {
        SetFirstException();
    }

    // TOCTOU race: We check exception state under lock, but call the callback without
    // holding the lock (to avoid deadlock). Another thread may set an exception between
    // the check and the callback invocation. This means onThreadExit may occasionally
    // execute despite an exception occurring. This is documented behavior - users requiring
    // strict "only on success" guarantees must implement their own exception tracking.
    // If the callback itself throws, it's caught and handled via SetFirstException() below.
    if (onThreadExit_) {
        std::unique_lock lock{exceptionMutex_};
        const bool okToCall{!exceptionOccurred_ && !stopSource_.stop_requested()};
        lock.unlock();

        if (okToCall) {
            try {
                onThreadExit_(threadIndex);
            } catch (...) {
                SetFirstException();
            }
        }
    }
}

template <typename ResultType>
ThreadPool<ResultType>::OptionalTask ThreadPool<ResultType>::PopTask(std::stop_token stopToken)
{
    std::unique_lock lock{poolMutex_};

    // Wait for work or shutdown
    pushed_.wait(lock, [this, &stopToken]() {
        if (stopToken.stop_requested() || stopSource_.stop_requested()) {
            return true;
        }

        const bool hasWork{!std::ranges::empty(head_)};
        const bool shuttingDown{!acceptingJobs_.load(std::memory_order_acquire)};

        if constexpr (IsFireAndForget) {
            // Fire-and-forget: wake on work or shutdown
            return hasWork || shuttingDown;
        } else {
            // Producer-consumer: wake on work AND tail space, or shutdown
            return (hasWork && (std::ranges::size(tail_) < maxQueueDepth_)) || shuttingDown;
        }
    });

    if (stopToken.stop_requested() || stopSource_.stop_requested() || std::ranges::empty(head_)) {
        return {};
    }

    // Extract highest priority task (min-heap via max-comparator)
    std::ranges::pop_heap(head_, TaskPriorityCompare{});
    PrioritizedTask& top{head_.back()};

    // For producer-consumer mode, enqueue future before executing task
    if constexpr (IsProducerConsumer) {
        tail_.push_back(top.Task.get_future());
        IncrementResultQueueDepth();
    }

    OptionalTask task{std::move(top.Task)};
    head_.pop_back();

    // Update metrics after dequeue
    UpdateMetricsOnQueueChange(std::ranges::size(head_));

    // Notify producer that space is available in input queue
    cvInputQueueNotFull_.notify_one();

    // Notify consumer that result is available (producer-consumer mode only)
    if constexpr (IsProducerConsumer) {
        cvResultQueueNotEmpty_.notify_one();
    }

    // Notify Finalize if head_ became empty (for fire-and-forget mode drain check)
    if (std::ranges::empty(head_)) {
        cvAllTasksComplete_.notify_all();
    }

    return task;
}

template <typename ResultType>
void ThreadPool<ResultType>::ExecuteTask(TaskType& task, ThreadIndex threadIndex)
{
    // Track active task count; RAII ensures decrement on all exit paths
    activeTasks_.fetch_add(1, std::memory_order_acq_rel);
    UpdateMetricsOnTaskStart();

    bool taskFailed{false};
    // If we captured taskFailed directly, we'll trigger -Wunused-lambda-capture
    const auto notifyOnExit = [&](void*) noexcept {
        const std::size_t prev{activeTasks_.fetch_sub(1, std::memory_order_release)};
        if constexpr (IsFireAndForget) {
            UpdateMetricsOnTaskComplete(taskFailed);
        } else {
            UpdateMetricsOnTaskComplete();
            (void)taskFailed;
        }
        if (prev == 1) {
            // Last task completed - wake threads waiting in Finalize()
            std::lock_guard lock{poolMutex_};
            cvAllTasksComplete_.notify_all();
        }
    };
    std::unique_ptr<void, decltype(notifyOnExit)> guard{this, notifyOnExit};

    if constexpr (IsFireAndForget) {
        // Fire-and-forget: get future to capture exceptions
        std::future<ResultType> future{task.get_future()};
        task(threadIndex);
        try {
            future.get();
        } catch (...) {
            taskFailed = true;
            throw;
        }
    } else {
        // Producer-consumer: future already in tail queue
        // packaged_task captures exceptions into the future, so we can't
        // detect them here. We'll detect them by checking the future's state.
        task(threadIndex);
    }
}

template <typename ResultType>
void ThreadPool<ResultType>::ProduceTask(TaskType&& task, std::int32_t priority)
{
    std::unique_lock lock{poolMutex_};

    if constexpr (IsProducerConsumer) {
        // Enforce consumer registration to prevent producer deadlock
        ThrowIfAborted();
        if (!acceptingJobs_.load(std::memory_order_acquire)) {
            throw std::runtime_error{"ThreadPool: not accepting jobs after finalization"};
        }
        cvConsumerRegistered_.wait(lock, [this] {
            return stopSource_.stop_requested() ||
                   !acceptingJobs_.load(std::memory_order_acquire) || consumerActive_;
        });
        ThrowIfAborted();
        if (!acceptingJobs_.load(std::memory_order_acquire)) {
            throw std::runtime_error{"ThreadPool: not accepting jobs after finalization"};
        }
        if (!consumerActive_) {
            throw std::runtime_error{
                "ThreadPool: consumer must call ConsumeWith() before submitting tasks"};
        }
    }

    // Block if work queue is full (back-pressure)
    // Predicate must not throw - check conditions and return true to wake up
    cvInputQueueNotFull_.wait(lock, [this] {
        // Return true if we should wake up to handle error/stop
        return stopSource_.stop_requested() || !acceptingJobs_.load(std::memory_order_acquire) ||
               (std::ranges::size(head_) < maxQueueDepth_);
    });

    // Check for abort/stop after wait returns (predicate must not throw)
    ThrowIfAborted();
    if (!acceptingJobs_.load(std::memory_order_acquire)) {
        throw std::runtime_error{"ThreadPool: not accepting jobs after finalization"};
    }

    // Enqueue with priority and sequence for stable priority queue
    const std::size_t sequence{nextSequence_.fetch_add(1, std::memory_order_relaxed)};
    head_.push_back(PrioritizedTask{priority, sequence, std::move(task)});
    std::ranges::push_heap(head_, TaskPriorityCompare{});

    // Update metrics after enqueue
    UpdateMetricsOnSubmit();
    UpdateMetricsOnQueueChange(std::ranges::size(head_));

    pushed_.notify_one();
}

template <typename ResultType>
void ThreadPool<ResultType>::SetFirstException(std::exception_ptr ptr) noexcept
{
    if (!ptr) {
        ptr = std::current_exception();
        if (!ptr) {
            // Defensive fallback: Should never happen in normal operation.
            // Only reachable if SetFirstException() is called incorrectly without
            // an active exception. Create a generic exception to avoid nullptr.
            ptr = std::make_exception_ptr(std::runtime_error{"ThreadPool: operation aborted"});
        }
    }

    {
        std::unique_lock lock{exceptionMutex_};
        if (exceptionOccurred_) {
            return;  // another thread beat us to it
        }

        exceptionOccurred_ = true;
        exc_ = std::move(ptr);
    }

    // Initiate fail-fast shutdown
    stopSource_.request_stop();

    // Drop queued tasks early to release captures on fail-fast
    {
        std::lock_guard lock{poolMutex_};
        head_.clear();
        UpdateMetricsOnQueueChange(0);
    }

    WakeAllWaiters();
}

template <typename ResultType>
void ThreadPool<ResultType>::ThrowIfAborted()
{
    if (stopSource_.stop_requested()) {
        std::lock_guard lock{exceptionMutex_};
        if (exc_) {
            std::rethrow_exception(exc_);
        }
    }
}

template <typename ResultType>
void ThreadPool<ResultType>::SignalNoMoreTasks() noexcept
{
    // Use release ordering: we're only publishing the "not accepting jobs" state
    if (acceptingJobs_.exchange(false, std::memory_order_release)) {
        {
            std::lock_guard lock{poolMutex_};
            pushed_.notify_all();
            cvConsumerRegistered_.notify_all();
        }
    }
}

template <typename ResultType>
void ThreadPool<ResultType>::JoinAllThreads(bool requestStop) noexcept
{
    if (requestStop) {
        stopSource_.request_stop();
        WakeAllWaiters();
    }

    // Copy thread IDs under lock to avoid race with IsWorkerThread().
    // We cannot move threads_ because IsWorkerThread() iterates over it.
    std::vector<std::thread::id> threadIds;
    try {
        std::lock_guard lock{poolMutex_};
        threadIds.reserve(std::ranges::size(threads_));
        for (const auto& thread : threads_) {
            threadIds.push_back(thread.get_id());
        }
    } catch (...) {
        // If lock fails, we can't safely access threads_, so we can't join them.
        // This should never happen in practice, but maintain noexcept guarantee.
        return;
    }

    // Join threads by ID (must re-acquire lock to access threads_)
    const std::thread::id self{std::this_thread::get_id()};
    for (const auto& threadId : threadIds) {
        if (threadId == self) {
            continue;  // never join self
        }
        // Find thread by ID and join it (release lock before joining)
        std::jthread* threadToJoin{nullptr};
        {
            std::lock_guard lock{poolMutex_};
            for (auto& thread : threads_) {
                if ((thread.get_id() == threadId) && thread.joinable()) {
                    threadToJoin = &thread;
                    break;
                }
            }
        }
        // Join outside lock to avoid deadlock
        if (threadToJoin && threadToJoin->joinable()) {
            threadToJoin->join();
        }
    }
}

template <typename ResultType>
template <typename Predicate>
void ThreadPool<ResultType>::WaitForCondition(Predicate&& pred)
{
    std::unique_lock lock{poolMutex_};
    cvAllTasksComplete_.wait(lock, std::forward<Predicate>(pred));
}

template <typename ResultType>
void ThreadPool<ResultType>::WakeAllWaiters() noexcept
{
    pushed_.notify_all();
    cvInputQueueNotFull_.notify_all();
    cvResultQueueNotEmpty_.notify_all();
    cvConsumerRegistered_.notify_all();
    cvAllTasksComplete_.notify_all();
}

template <typename F, typename PoolType>
void Dispatch(PoolType&& pool, F&& callback, std::int32_t numEntries,
              std::optional<std::int32_t> priority = std::nullopt)
{
    using T = std::remove_cvref_t<PoolType>;
    ThreadPool<>* tp{nullptr};

    if constexpr (std::same_as<T, ThreadPool<>>) {
        tp = std::addressof(pool);
    } else if constexpr (std::is_pointer_v<T>) {
        tp = pool;
    } else if constexpr (std::same_as<T, std::reference_wrapper<ThreadPool<>>>) {
        tp = std::addressof(pool.get());
    } else if constexpr (std::same_as<T, std::optional<ThreadPool<>>>) {
        tp = pool ? &*pool : nullptr;
    } else if constexpr (std::same_as<T, std::optional<ThreadPool<>*>>) {
        tp = pool ? *pool : nullptr;
    } else if constexpr (std::same_as<T, std::optional<std::reference_wrapper<ThreadPool<>>>>) {
        tp = pool ? &pool->get() : nullptr;
    } else if constexpr (requires {
                             { pool.get() } -> std::same_as<ThreadPool<>*>;
                         }) {
        tp = pool.get();
    }

    if (tp) {
        if (priority) {
            tp->DispatchWithPriority(numEntries, *priority, std::forward<F>(callback));
        } else {
            tp->Dispatch(numEntries, std::forward<F>(callback));
        }
    } else {
        // No pool: execute serially. Apply same validation as member function.
        if (numEntries < 0) {
            throw std::invalid_argument{"Dispatch: numEntries must be >= 0"};
        }
        for (std::int32_t i{0}; i < numEntries; ++i) {
            callback(i);
        }
    }
}

}  // namespace Parallel
}  // namespace PacBio

#endif  // PBCOPPER_PARALLEL_THREADPOOL_H
