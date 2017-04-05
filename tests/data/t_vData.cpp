#include <catch.hpp>
#include "data/vData.hpp"

using namespace Anaquin;

TEST_CASE("ReadVFile_Example")
{
    const auto r = readVFile(Reader("data/VarQuin/AVA009_v001.vcf"));
    
    REQUIRE(r.countInd()    == 108);
    REQUIRE(r.countIndSyn() == 108);
    REQUIRE(r.countIndGen() == 0);
    REQUIRE(r.countSNP()    == 137);
    REQUIRE(r.countSNPSyn() == 137);
    REQUIRE(r.countSNPGen() == 0);    
}
