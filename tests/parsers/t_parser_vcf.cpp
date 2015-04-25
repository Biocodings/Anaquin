#include <catch.hpp>
#include "parser_vcf.hpp"

using namespace Spike;

TEST_CASE("DNA_Variant")
{
    ParserVCF::parse("tests/data/dna_sims/DNA.flat.chrT.vcf", [&](const VCFVariant &v)
    {
        // Empty Implementation        
    });

    //ASSERT_EQ(0, i);
    //ASSERT_EQ(245, j);
}
