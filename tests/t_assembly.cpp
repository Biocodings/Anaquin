#include "gtest/gtest.h"
#include "assembly.hpp"

TEST(Generated, AssemblyTest)
{
    const auto r = Assembly::analyze("/Users/tedwong/Sources/ABCD/transcripts/transcripts.gtf");

    ASSERT_EQ(355, r.exon.tp);
    ASSERT_EQ(21, r.exon.fp);
}