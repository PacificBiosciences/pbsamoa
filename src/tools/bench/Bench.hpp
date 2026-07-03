#ifndef PBSAMOA_TOOLS_BENCH_HPP
#define PBSAMOA_TOOLS_BENCH_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace Bench {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace Bench
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_BENCH_HPP
