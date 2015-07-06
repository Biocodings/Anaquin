#include <catch.hpp>
#include "fusion/f_fusion.hpp"

using namespace Spike;

TEST_CASE("FFusion_Test")
{
    const auto stats = FFusion::analyze("tests/data/fusion/fusions.out");
    
    REQUIRE(stats.p.m.sn() == Approx(0.125));
    REQUIRE(stats.p.m.sp() == Approx(0.3157894737));
    REQUIRE(stats.p.m.nr == 48);
}