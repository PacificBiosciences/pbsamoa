#ifndef PBSAMOA_TOOLS_ZMI_INDEX_HPP
#define PBSAMOA_TOOLS_ZMI_INDEX_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace ZmiIndex {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace ZmiIndex
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_ZMI_INDEX_HPP
