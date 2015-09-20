#include <catch.hpp>
#include "unit/test.hpp"
#include "fusion/f_discover.hpp"

using namespace Anaquin;

TEST_CASE("FDiscover_F_1001")
{
    Test::fusionA();
    
    const auto r = Test::test("-t FusionDiscover -rfus tests/data/F_1000/FUSStandard_1.0.ref -m tests/data/F_1000/FUSMixture_1.0.csv -soft tophat -fuzzy 100 -uout tests/data/F_1001/fusions.out");
    
    REQUIRE(r.status == 0);
    
    Test::fusionA();

    const auto stats = FDiscover::report("tests/data/F_1001/fusions.out", FDiscover::Options(FDiscover::Software::TopHat, 100));

    REQUIRE(stats.m.sn() == Approx(0.9166666667));
}

TEST_CASE("FDiscover_F_1000")
{
    Test::fusionA();

    const auto r = Test::test("-t FusionDiscover -rfus tests/data/F_1000/FUSStandard_1.0.ref -m tests/data/F_1000/FUSMixture_1.0.csv -soft star -uout tests/data/F_1000/star-fusion.fusion_candidates.txt");
    
    REQUIRE(r.status == 0);
    
    Test::fusionA();

    const auto stats = FDiscover::report("tests/data/F_1000/star-fusion.fusion_candidates.txt",
                                                FDiscover::Options(FDiscover::Software::Star));
    const auto lm = stats.linear();

    REQUIRE(lm.r  == Approx(0.9729157505));
    REQUIRE(lm.m  == Approx(0.9771885928));
    REQUIRE(lm.r2 == Approx(0.9465650576));
}