#include "../../src/ZmwUtils.hpp"

#include <pbsamoa/core/Tags.hpp>

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

TEST(ZmwUtils, MovieZmwPrefixTrailingSecondFieldSlash)
{
    EXPECT_EQ(MovieZmwPrefix("movie/42/"), "movie/42");
}

TEST(ZmwUtils, MovieZmwPrefixEmptySecondField)
{
    EXPECT_EQ(MovieZmwPrefix("movie//rest"), "movie/");
}

TEST(ZmwUtils, MovieZmwPrefixDifferentMoviesSameZmw)
{
    EXPECT_NE(MovieZmwPrefix("movieA/42/ccs"), MovieZmwPrefix("movieB/42/ccs"));
}

TEST(ZmwUtils, MovieZmwPrefixSameMovieDifferentZmw)
{
    EXPECT_NE(MovieZmwPrefix("movie/42/ccs"), MovieZmwPrefix("movie/99/ccs"));
}

// --- ReadGroupBaseId / ParseReadGroupId ---

TEST(ZmwUtils, ReadGroupBaseIdWithoutBarcodeSuffix)
{
    EXPECT_EQ(ReadGroupBaseId("12345678"), "12345678");
}

TEST(ZmwUtils, ReadGroupBaseIdWithBarcodeSuffix)
{
    EXPECT_EQ(ReadGroupBaseId("12345678/0--0"), "12345678");
}

TEST(ZmwUtils, ParseReadGroupIdHex)
{
    EXPECT_EQ(ParseReadGroupId("0000002a"), 42);
    EXPECT_EQ(ParseReadGroupId("0000002a/0--0"), 42);
}

TEST(ZmwUtils, ParseReadGroupIdInvalidReturnsZero)
{
    EXPECT_EQ(ParseReadGroupId("rg1"), 0);
    EXPECT_EQ(ParseReadGroupId(""), 0);
}

// --- ParseZmwIdentity ---

TEST(ZmwUtils, ParseZmwIdentityNormal)
{
    TagMap tags;
    tags.Set(RG_TAG, std::string{"0000002a/0--0"});
    EXPECT_EQ(ParseZmwIdentity("movie/42/0_100", tags), (ZmwIdentity{42, 42}));
}

TEST(ZmwUtils, ParseZmwIdentityMissingRgFallsBackToZero)
{
    const TagMap tags;
    EXPECT_EQ(ParseZmwIdentity("movie/42/0_100", tags), (ZmwIdentity{0, 42}));
}

TEST(ZmwUtils, ParseZmwIdentityInvalidRgTextFallsBackToZero)
{
    TagMap tags;
    tags.Set(RG_TAG, std::string{"not-hex"});
    EXPECT_EQ(ParseZmwIdentity("movie/42/0_100", tags), (ZmwIdentity{0, 42}));
}

}  // namespace Samoa
}  // namespace PacBio
