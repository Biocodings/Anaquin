#include "gtest/gtest.h"
#include "AssemblyAnalyst.hpp"

TEST(Generated, AssemblyCompareTest)
{
    const auto r = AssemblyAnalyst::analyze("/Users/tedwong/Sources/ABCD/transcripts/transcripts.gtf");

    ASSERT_EQ(355, r.exon.tp);
    ASSERT_EQ(21, r.exon.fp);
}