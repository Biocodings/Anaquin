#include <catch.hpp>
#include "tools/gtf_data.hpp"

using namespace Anaquin;

TEST_CASE("GTF_Summary_1")
{
    const auto r = gtfData(Reader("data/RnaQuin/ATR001.v032.gtf"));
    const auto i = r.gIntervals(ChrT);
    
    REQUIRE(i.size() == 78);
    REQUIRE(i.overlap(Locus(6955490, 6955495)));
    REQUIRE(i.contains(Locus(6955490, 6955495)));
    REQUIRE(!i.contains(Locus(6955480, 6955485)));
    REQUIRE(!i.overlap(Locus(6955480, 6955485)));
    
    REQUIRE(r.countGene()     == 78);
    REQUIRE(r.countGeneSyn()  == 78);
    REQUIRE(r.countTrans()    == 164);
    REQUIRE(r.countTransSyn() == 164);
    REQUIRE(r.countExon()     == 1192);
    REQUIRE(r.countExonSyn()  == 1192);
    REQUIRE(r.countUExon()    == 869);
    REQUIRE(r.countUExonSyn() == 869);
    REQUIRE(r.countIntr()     == 1028);
    REQUIRE(r.countIntrSyn()  == 1028);
    REQUIRE(r.countUIntr()    == 754);
    REQUIRE(r.countUIntrSyn() == 754);

    REQUIRE(r.il.at(ChrT).at(Locus(6955730, 6960383)) == 1);
    REQUIRE(r.il.at(ChrT).at(Locus(2227518, 2235700)) == 3);
}

TEST_CASE("GTF_Summary_2")
{
    const auto r = gtfData(Reader("data/tests/RnaQuin/combined.gtf"));
    
    REQUIRE(r.countGene()     == 958);
    REQUIRE(r.countGeneSyn()  == 78);
    REQUIRE(r.countTrans()    == 2599);
    REQUIRE(r.countTransSyn() == 164);
    REQUIRE(r.countExon()     == 15110);
    REQUIRE(r.countExonSyn()  == 1192);
    REQUIRE(r.countIntr()     == 12511);
    REQUIRE(r.countIntrSyn()  == 1028);
}