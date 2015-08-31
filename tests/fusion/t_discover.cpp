#include <catch.hpp>
#include "unit/test.hpp"
#include "fusion/f_discover.hpp"

using namespace Anaquin;

TEST_CASE("FDiscover_F_1000_Command")
{
    const auto r = Test::test("-t FusionDiscover -rfus tests/data/F_1000/FUSStandard_1.0.ref -m tests/data/F_1000/FUSMixture_1.0.csv -soft star -uout tests/data/F_1000/star-fusion.fusion_candidates.txt");
    
    REQUIRE(r.status == 0);
}

TEST_CASE("FDiscover_F_1000")
{
    Standard::instance(true);
    Standard::instance().f_ref(Reader("tests/data/F_1000/FUSStandard_1.0.ref"));
    Standard::instance().f_mix(Reader("tests/data/F_1000/FUSMixture_1.0.csv"));
    
    const auto stats = FDiscover::analyze("tests/data/F_1000/star-fusion.fusion_candidates.txt",
                                                FDiscover::Options(Software::Star));
    const auto lm = stats.cov.linear();

    REQUIRE(lm.r  == Approx(0.9729157505));
    REQUIRE(lm.m  == Approx(0.9771885928));
    REQUIRE(lm.r2 == Approx(0.9434218257));
}

//TEST_CASE("FDiscover_Simulated")
//{
//    const auto r = Test::test("-t FusionDiscover -rfus data/fusion/FUSBreak_1.0.ref -m data/fusion/FUSMixture_3.0.csv -uout tests/data/fusion/simulated/fusions.out");
//    
//    REQUIRE(r.status == 0);
//    REQUIRE(r.output.find("Fusion Analysis") != std::string::npos);
//
//    // Apply default resources
//    Test::fusion();
//
//    const auto stats = FDiscover::analyze("tests/data/fusion/simulated/fusions.out");
//    
//    // The linear model associated with the expression
//    const auto lm = stats.linear();
//    
//    REQUIRE(lm.r  == Approx(0.9827924165));
//    REQUIRE(lm.r2 == Approx(0.96387393));
//    REQUIRE(lm.m  == Approx(0.8617847724));
//    REQUIRE(lm.c  == Approx(3.2714500279));
//}
//
//TEST_CASE("FDiscover_10K")
//{
//    const auto r = Test::test("-t FusionDiscover -rfus data/fusion/FUSBreak_1.0.ref -m data/fusion/FUSMixture_3.0.csv -uout tests/data/fusion/10K/fusions.out");
//
//    REQUIRE(r.status == 0);
//    REQUIRE(r.output.find("Fusion Analysis") != std::string::npos);
//}