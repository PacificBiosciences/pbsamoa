#ifndef PBSAMOA_TOOLS_BAIQUERY_HPP
#define PBSAMOA_TOOLS_BAIQUERY_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace BaiQuery {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace BaiQuery
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_BAIQUERY_HPP
