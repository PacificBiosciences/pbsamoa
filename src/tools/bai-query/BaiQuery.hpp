#ifndef PBSAMOA_TOOLS_BAIQUERY_HPP
#define PBSAMOA_TOOLS_BAIQUERY_HPP

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Results.h>

namespace PacBio {
namespace Samoa {
namespace BaiQuery {

CLI_v2::Interface CreateInterface();
int Runner(const CLI_v2::Results& results);

}  // namespace BaiQuery
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_BAIQUERY_HPP
