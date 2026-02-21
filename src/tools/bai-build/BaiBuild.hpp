#ifndef PBSAMOA_TOOLS_BAIBUILD_HPP
#define PBSAMOA_TOOLS_BAIBUILD_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace BaiBuild {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace BaiBuild
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_BAIBUILD_HPP
