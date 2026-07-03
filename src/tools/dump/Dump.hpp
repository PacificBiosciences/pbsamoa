#ifndef PBSAMOA_TOOLS_DUMP_HPP
#define PBSAMOA_TOOLS_DUMP_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace Dump {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace Dump
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_DUMP_HPP
