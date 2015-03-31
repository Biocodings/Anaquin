#include <catch.hpp>
#include "parser_vcf.hpp"

TEST_CASE("DNA_Variant")
{
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/variant.ChrT51.vcf", [&](const VCFVariant &v)
    {
    });

    //ASSERT_EQ(0, i);
    //ASSERT_EQ(245, j);
}

TEST_CASE("DNA_Hetero")
{
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/hetero.ChrT51.vcf", [&](const VCFVariant &v)
                     {
                     });
}

TEST_CASE("DNA_Ref")
{
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/ref.ChrT51.vcf", [&](const VCFVariant &v)
                     {
                     });
}