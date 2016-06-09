#include <catch.hpp>
#include "RnaQuin/r_viewer.hpp"

using namespace Anaquin;

#ifdef REGRESSION

TEST_CASE("TViewer_T_1001")
{
    const auto r1 = Test::test("-t TransIGV -o /data/regress/latest/T_1001 -ubam /share/Projects/TransQuin/K_RMXA/A1/accepted_hits.bam");
    const auto r2 = Test::test("-t TransIGV -o /data/regress/latest/T_1001 -ubam /share/Projects/TransQuin/K_RMXA/A2/accepted_hits.bam");
    const auto r3 = Test::test("-t TransIGV -o /data/regress/latest/T_1001 -ubam /share/Projects/TransQuin/K_RMXA/A3/accepted_hits.bam");
    
    REQUIRE(r1.status == 0);
    REQUIRE(r2.status == 0);
    REQUIRE(r3.status == 0);
}

#endif
