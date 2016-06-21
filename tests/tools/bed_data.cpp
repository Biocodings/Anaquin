#include <catch.hpp>
#include "tools/bed_data.hpp"

using namespace Anaquin;

TEST_CASE("BED_Synthetic")
{
    const auto r = bedData(Reader("data/VarQuin/AVA017.v032.bed"));
    
    REQUIRE(r.countGene()    == 36);
    REQUIRE(r.countGeneSyn() == 36);
    REQUIRE(r.countGeneGen() == 0);
    
    const auto i = r.inters();
    
    REQUIRE(i.size() == 1);
    
    REQUIRE(i.at(ChrT).exact(Locus(373692, 374677)));
    REQUIRE(i.at(ChrT).contains(Locus(373692, 374677)));
    REQUIRE(i.at(ChrT).overlap(Locus(373692, 374677)));
    
    REQUIRE(!i.at(ChrT).exact(Locus(373691, 374677)));
    REQUIRE(!i.at(ChrT).contains(Locus(373691, 374677)));
    REQUIRE(i.at(ChrT).overlap(Locus(373691, 374677)));
}