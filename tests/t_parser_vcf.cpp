#include <catch.hpp>
#include "parser_vcf.hpp"

TEST_CASE("DNA_Variant")
{
    ParserVCF::parse("data/DNA/variant.ChrT51.vcf", [&](const VCFVariant &v)
    {
    });

    //ASSERT_EQ(0, i);
    //ASSERT_EQ(245, j);
}
