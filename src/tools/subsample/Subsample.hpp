#ifndef PBSAMOA_TOOLS_SUBSAMPLE_HPP
#define PBSAMOA_TOOLS_SUBSAMPLE_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace Subsample {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace Subsample
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_SUBSAMPLE_HPP
