#include <catch.hpp>
#include "dna/d_align.hpp"

using namespace Spike;

TEST_CASE("DAlign_Simulations")
{
    const auto r = DAlign::analyze("tests/data/dna_sims/aligned.sam");

    REQUIRE(r.me.sp() == Approx(0.999984782));
    REQUIRE(r.me.sn() == Approx(0.9997261462));
    REQUIRE(r.me.nq == 262846);
    REQUIRE(r.me.nr == 262914);
    //REQUIRE(r.se.id == "R_5_1");
    //REQUIRE(r.se.counts == 2);
    //REQUIRE(r.se.abund == 9.765625);
}