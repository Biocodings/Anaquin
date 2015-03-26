#include "gtest/gtest.h"
#include "parser_gtf.hpp"

TEST(Options, ParserGTFTest)
{
    ParserGTF::parse("/Users/tedwong/Sources/ABCD/standards/RNAstandards.gtf", [&](const Feature &f, ParserProgress &p)
    {
        ASSERT_EQ(2, f.options.size());
        ASSERT_TRUE(f.options.count("gene_id"));
        ASSERT_TRUE(f.options.count("transcript_id"));
    });
}