#include <catch.hpp>
#include "parser_vcf.hpp"

TEST_CASE("DNA_Variant")
{
    ParserVCF::parse("tests/data/simulation/DNA.flat.chrT.vcf", [&](const VCFVariant &v)
    {
        
        
    });

    //ASSERT_EQ(0, i);
    //ASSERT_EQ(245, j);
}
