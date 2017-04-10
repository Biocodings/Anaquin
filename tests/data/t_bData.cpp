#include <catch.hpp>
#include "data/bData.hpp"
#include "RnaQuin/RnaQuin.hpp"

using namespace Anaquin;

TEST_CASE("ReadRegions_Example")
{
    const auto r = readRegions(Reader("data/VarQuin/AVA017_v001.bed"));
    
    auto allTrue = [&](const ChrID &)
    {
        return true;
    };
    
    REQUIRE(r.nGene()    == 36);
    REQUIRE(r.nGeneSyn(allTrue) == 36);
    REQUIRE(r.nGeneGen(allTrue) == 0);

    const auto i = r.inters();
    
    REQUIRE(i.size() == 1);
    
    REQUIRE(i.at(ChrIS).exact(Locus(373692, 374677)));
    REQUIRE(i.at(ChrIS).contains(Locus(373692, 374677)));
    REQUIRE(i.at(ChrIS).overlap(Locus(373692, 374677)));
    
    REQUIRE(!i.at(ChrIS).exact(Locus(373691, 374677)));
    REQUIRE(!i.at(ChrIS).contains(Locus(373691, 374677)));
    REQUIRE(i.at(ChrIS).overlap(Locus(373691, 374677)));
}
