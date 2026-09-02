#ifndef PBSAMOA_TOOLS_UNDO_HPP
#define PBSAMOA_TOOLS_UNDO_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace Undo {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace Undo
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_UNDO_HPP
