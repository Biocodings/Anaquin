#include <catch.hpp>
#include "tools/bed_data.hpp"

using namespace Anaquin;

TEST_CASE("BED_Synthetic")
{
    const auto r = bedData(Reader("data/VarQuin/AVA017_v001.bed"));
    
    REQUIRE(r.nGene()    == 36);
    REQUIRE(r.nGeneSyn() == 36);
    REQUIRE(r.nGeneGen() == 0);
    
    const auto i = r.inters();
    
    REQUIRE(i.size() == 1);
    
    REQUIRE(i.at(ChrIS).exact(Locus(373692, 374677)));
    REQUIRE(i.at(ChrIS).contains(Locus(373692, 374677)));
    REQUIRE(i.at(ChrIS).overlap(Locus(373692, 374677)));
    
    REQUIRE(!i.at(ChrIS).exact(Locus(373691, 374677)));
    REQUIRE(!i.at(ChrIS).contains(Locus(373691, 374677)));
    REQUIRE(i.at(ChrIS).overlap(Locus(373691, 374677)));
}
