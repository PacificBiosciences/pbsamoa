#ifndef PBSAMOA_TOOLS_BENCH_HPP
#define PBSAMOA_TOOLS_BENCH_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace Bench {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace Bench
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_BENCH_HPP
