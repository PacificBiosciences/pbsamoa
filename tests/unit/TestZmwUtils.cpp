#include "../../src/ZmwUtils.hpp"

#include <gtest/gtest.h>

#include <string_view>

namespace PacBio {
namespace Samoa {

// --- ParseZmwFromName ---

TEST(ZmwUtils, ParseZmwFromNameNormal)
{
    EXPECT_EQ(ParseZmwFromName("movie/42/0_100"), 42);
    EXPECT_EQ(ParseZmwFromName("m12345_20250101/123456/ccs"), 123456);
}

TEST(ZmwUtils, ParseZmwFromNameNoSlash) { EXPECT_EQ(ParseZmwFromName("noslash"), 0); }

TEST(ZmwUtils, ParseZmwFromNameOneSlash) { EXPECT_EQ(ParseZmwFromName("movie/42"), 42); }

TEST(ZmwUtils, ParseZmwFromNameEmptyZmw) { EXPECT_EQ(ParseZmwFromName("movie//rest"), 0); }

TEST(ZmwUtils, ParseZmwFromNameNonNumericZmw) { EXPECT_EQ(ParseZmwFromName("movie/abc/rest"), 0); }

TEST(ZmwUtils, ParseZmwFromNameEmpty) { EXPECT_EQ(ParseZmwFromName(""), 0); }

// --- MovieZmwPrefix ---

TEST(ZmwUtils, MovieZmwPrefixNormal)
{
    EXPECT_EQ(MovieZmwPrefix("movie/42/0_100"), "movie/42");
    EXPECT_EQ(MovieZmwPrefix("m12345/999/ccs"), "m12345/999");
}

TEST(ZmwUtils, MovieZmwPrefixNoSlash) { EXPECT_EQ(MovieZmwPrefix("noslash"), "noslash"); }

TEST(ZmwUtils, MovieZmwPrefixOneSlash) { EXPECT_EQ(MovieZmwPrefix("movie/42"), "movie/42"); }

TEST(ZmwUtils, MovieZmwPrefixDifferentMoviesSameZmw)
{
    EXPECT_NE(MovieZmwPrefix("movieA/42/ccs"), MovieZmwPrefix("movieB/42/ccs"));
}

TEST(ZmwUtils, MovieZmwPrefixSameMovieDifferentZmw)
{
    EXPECT_NE(MovieZmwPrefix("movie/42/ccs"), MovieZmwPrefix("movie/99/ccs"));
}

}  // namespace Samoa
}  // namespace PacBio
