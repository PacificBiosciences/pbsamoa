#ifndef PBSAMOA_TOOLS_ZMIQUERY_HPP
#define PBSAMOA_TOOLS_ZMIQUERY_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace ZmiQuery {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace ZmiQuery
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_ZMIQUERY_HPP
