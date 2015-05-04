#include <catch.hpp>
#include "parsers/parser_vcf.hpp"

using namespace Spike;

TEST_CASE("DNA_Variant")
{
    ParserVCF::parse("tests/data/dna_sims/DNA.flat.chrT.vcf", [&](const VCFVariant &v, const ParserProgress &)
    {
        // Empty Implementation        
    });

    //ASSERT_EQ(0, i);
    //ASSERT_EQ(245, j);
}
