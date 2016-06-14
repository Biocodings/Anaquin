#include <catch.hpp>
#include "tools/summary.hpp"

using namespace Anaquin;

TEST_CASE("GTF_Summary_1")
{
    const auto r = summaryGTF(Reader("data/RnaQuin/ATR001.v032.gtf"));
    
    REQUIRE(r.countGenes() == 78);
    REQUIRE(r.countTrans() == 164);
    REQUIRE(r.countExons() == 1192);
    REQUIRE(r.countIntrs() == 1028);
}

TEST_CASE("GTF_Summary_2")
{
    const auto r = summaryGTF(Reader("data/tests/RnaQuin/combined.gtf"));
    
    REQUIRE(r.countGenes() == 958);
    REQUIRE(r.countTrans() == 2599);
    REQUIRE(r.countExons() == 15110);
    REQUIRE(r.countIntrs() == 12511);
}