#ifndef PBSAMOA_TOOLS_MERGE_HPP
#define PBSAMOA_TOOLS_MERGE_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace MergeTool {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace MergeTool
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_MERGE_HPP
