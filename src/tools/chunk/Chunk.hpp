#ifndef PBSAMOA_TOOLS_CHUNK_HPP
#define PBSAMOA_TOOLS_CHUNK_HPP

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Samoa {
namespace ChunkTool {

PacBio::CLI_v2::Interface CreateInterface();
int Runner(const PacBio::CLI_v2::Results& results);

}  // namespace ChunkTool
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_CHUNK_HPP
