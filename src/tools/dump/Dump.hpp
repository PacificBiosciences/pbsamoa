#ifndef PBSAMOA_TOOLS_DUMP_HPP
#define PBSAMOA_TOOLS_DUMP_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace Dump {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace Dump
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_DUMP_HPP
