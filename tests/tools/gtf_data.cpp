#include <catch.hpp>
#include "tools/gtf_data.hpp"

using namespace Anaquin;

TEST_CASE("GTF_Summary_1")
{
    const auto r = gtfData(Reader("data/RnaQuin/ATR001.v032.gtf"));
    
    REQUIRE(r.countGene()  == 78);
    REQUIRE(r.countTrans() == 164);
    REQUIRE(r.countExon()  == 1192);
    REQUIRE(r.countIntr()  == 1028);
}

TEST_CASE("GTF_Summary_2")
{
    const auto r = gtfData(Reader("data/tests/RnaQuin/combined.gtf"));
    
    REQUIRE(r.countGene()  == 958);
    REQUIRE(r.countTrans() == 2599);
    REQUIRE(r.countExon()  == 15110);
    REQUIRE(r.countIntr()  == 12511);
}