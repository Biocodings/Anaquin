#include "gtest/gtest.h"
#include "ParserGTF.hpp"

using namespace std;

TEST(Options, ParserGTFTest)
{
    ParserGTF::parse("/Users/user1/Sources/ABCD/standards/RNAstandards.gtf", [&](const Feature &f, ParserProgress &p)
    {
        ASSERT_EQ(2, f.options.size());
        ASSERT_TRUE(f.options.count("gene_id"));
        ASSERT_TRUE(f.options.count("transcript_id"));
    });
}