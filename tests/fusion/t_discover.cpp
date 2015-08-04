#include <catch.hpp>
#include "unit/test.hpp"
#include "fusion/f_discover.hpp"

using namespace Anaquin;

TEST_CASE("FDiscover_100K_Star")
{
    FDiscover::Options o;

    // The experiment was run on FusionStar
    o.soft = Software::Star;

    const auto stats = FDiscover::analyze("tests/data/fusion/100K/star-fusion.fusion_candidates.txt", o);
}

TEST_CASE("FDiscover_Simulated")
{
    const auto stats = FDiscover::analyze("tests/data/fusion/simulated/fusions.out");
    
    // The linear model associated with the expression
    const auto lm = stats.linear();
    
    REQUIRE(lm.r  == Approx(0.9489461887));
    REQUIRE(lm.r2 == Approx(0.8960126545));
    REQUIRE(lm.m  == Approx(0.945210027));
    REQUIRE(lm.c  == Approx(2.5271092606));

    const auto cmd = "-t FusionDiscover -rfus data/fusion/FUS.v1.ref -m data/fusion/FUS.v3.csv -uout tests/data/fusion/simulated/fusions.out";
    
    const auto r = CommandTest::test(cmd);

    REQUIRE(r.status == 0);
    REQUIRE(r.output.find("Fusion Analysis") != std::string::npos);
}

TEST_CASE("FDiscover_10K")
{
    const auto cmd = "-t FusionDiscover -rfus tests/data/fusion/10K/FUS.v1.ref -m tests/data/fusion/10K/FUSE_mixtures_v3.csv -uout tests/data/fusion/10K/fusions.out";

    const auto r = CommandTest::test(cmd);
    
    REQUIRE(r.status == 0);
    REQUIRE(r.output.find("Fusion Analysis") != std::string::npos);
}