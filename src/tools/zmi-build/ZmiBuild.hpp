#ifndef PBSAMOA_TOOLS_ZMIBUILD_HPP
#define PBSAMOA_TOOLS_ZMIBUILD_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace ZmiBuild {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace ZmiBuild
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_ZMIBUILD_HPP
