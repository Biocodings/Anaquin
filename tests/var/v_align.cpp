#include <catch.hpp>
#include "var/v_align.hpp"

using namespace Anaquin;

TEST_CASE("VAlign_Simulation")
{
    //const auto r = VAlign::analyze("tests/data/var/aligned.sam");
    const auto r = VAlign::analyze("ABCD.bam");

    REQUIRE(r.p.m.sp() == Approx(0.999984782));
    REQUIRE(r.p.m.sn() == Approx(0.9998288238));
    REQUIRE(r.p.m.nq == 262846);
    REQUIRE(r.p.m.nr == 262887);
    REQUIRE(r.p.s.id == "D_3_7_R");
    REQUIRE(r.p.s.counts == 2);
    REQUIRE(r.p.s.abund == Approx(1.005944252));
}