#include <catch.hpp>
#include "tools/summary.hpp"

using namespace Anaquin;

TEST_CASE("GTF_Summary_1")
{
    const auto r = summaryGTF(Reader("data/RnaQuin/ATR001.v032.gtf"));

    REQUIRE(r.countGenes() == 78);
    REQUIRE(r.countTrans() == 164);
}