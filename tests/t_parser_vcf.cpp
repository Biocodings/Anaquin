#include "gtest/gtest.h"
#include "parser_vcf.hpp"

TEST(DNA_Variant, ParserVCFTest)
{
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/variant.ChrT51.vcf", [&](const VCFVariant &v)
    {
    });

    //ASSERT_EQ(0, i);
    //ASSERT_EQ(245, j);
}

TEST(DNA_Hetero, ParserVCFTest)
{
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/hetero.ChrT51.vcf", [&](const VCFVariant &v)
                     {
                     });
}

TEST(DNA_Ref, ParserVCFTest)
{
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/ref.ChrT51.vcf", [&](const VCFVariant &v)
                     {
                     });
}