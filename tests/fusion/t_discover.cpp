#include <catch.hpp>
#include "unit/test.hpp"
#include "fusion/f_discover.hpp"

using namespace Anaquin;

TEST_CASE("FDiscover_F_1001_Invalid")
{
    Test::fusionA();
    
    const auto r = Test::test("-t FusionDiscover -rfus data/fusion/AFU004.v032.ref -soft tohat -fuzzy 100 -uout tests/data/F_1001/fusions.out");
    
    REQUIRE(r.status == 1);
}

TEST_CASE("FDiscover_F_1001")
{
    Test::fusionA();
    
    const auto r = Test::test("-t FusionDiscover -rfus data/fusion/AFU004.v032.ref -soft tophat -fuzzy 100 -uout tests/data/F_1001/fusions.out");
    
    REQUIRE(r.status == 0);
    
    Test::fusionA();

    const auto stats = FDiscover::report("tests/data/F_1001/fusions.out", FDiscover::Options(FDiscover::Software::TopHat, 100));

    REQUIRE(stats.m.sn() == Approx(0.9166666667));
}