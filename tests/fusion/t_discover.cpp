#include <catch.hpp>
#include "unit/test.hpp"
#include "fusion/f_discover.hpp"

using namespace Anaquin;

TEST_CASE("FDiscover_F_1001")
{
    Test::fusion();
    
    const auto r = Test::test("-t FusionDiscover -rfus tests/data/F_1000/FUSStandard_1.0.ref -m tests/data/F_1000/FUSMixture_1.0.csv -soft tophat -fuzzy 100 -uout tests/data/F_1001/fusions.out");
    
    REQUIRE(r.status == 0);
    
    Test::fusion();

    const auto stats = FDiscover::analyze("tests/data/F_1001/fusions.out", FDiscover::Options(Software::TopHat, 100));

    REQUIRE(stats.m.sn() == Approx(0.9166666667));
}

TEST_CASE("FDiscover_F_1000")
{
    Test::fusion();

    const auto r = Test::test("-t FusionDiscover -rfus tests/data/F_1000/FUSStandard_1.0.ref -m tests/data/F_1000/FUSMixture_1.0.csv -soft star -uout tests/data/F_1000/star-fusion.fusion_candidates.txt");
    
    REQUIRE(r.status == 0);
    
    Test::fusion();

    const auto stats = FDiscover::analyze("tests/data/F_1000/star-fusion.fusion_candidates.txt",
                                                FDiscover::Options(Software::Star));
    const auto lm = stats.cov.linear();

    REQUIRE(lm.r  == Approx(0.9729157505));
    REQUIRE(lm.m  == Approx(0.9771885928));
    REQUIRE(lm.r2 == Approx(0.9434218257));
}