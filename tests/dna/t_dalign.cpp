#include <catch.hpp>
#include "dna/d_align.hpp"

using namespace Spike;

TEST_CASE("DAlign_Simulations")
{
    const auto r = DAlign::analyze("/Users/tedwong/Sources/QA/flat.sam");

    REQUIRE(r.p.m.sp() == Approx(0.999984782));
    REQUIRE(r.p.m.sn() == Approx(0.9997261462));
    REQUIRE(r.p.m.nq == 262846);
    REQUIRE(r.p.m.nr == 262914);
    //REQUIRE(r.se.id == "R_5_1");
    //REQUIRE(r.se.counts == 2);
    //REQUIRE(r.se.abund == 9.765625);
}