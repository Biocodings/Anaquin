#include <catch.hpp>
#include "unit/test.hpp"

using namespace Anaquin;

TEST_CASE("Anaquin_Version")
{
    const auto r = Test::test("-v");
    
    REQUIRE(r.status == 0);
    REQUIRE(r.output == "Anaquin v1.1.1\n");
}

TEST_CASE("Anaquin_InvalidMixture")
{
    const auto r = Test::test("-t VarAllele -soft gatk -rvcf data/VarQuin/AVA009.v032.vcf -m variant.vcf -ufiles variant.vcf");
    
    
    
    
}

TEST_CASE("Anaquin_EmptyMixture")
{
    
}