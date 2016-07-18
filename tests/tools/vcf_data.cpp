#include <catch.hpp>
#include "tools/vcf_data.hpp"

using namespace Anaquin;

TEST_CASE("VCF_Synthetic")
{
    const auto r = vcfData(Reader("data/VarQuin/AVA009_v001.vcf"));
    
    REQUIRE(r.countInd()    == 108);
    REQUIRE(r.countIndSyn() == 108);
    REQUIRE(r.countIndGen() == 0);
    REQUIRE(r.countSNP()    == 137);
    REQUIRE(r.countSNPSyn() == 137);
    REQUIRE(r.countSNPGen() == 0);    
}