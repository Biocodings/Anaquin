#include <catch.hpp>
#include "fus/f_discover.hpp"

using namespace Spike;

TEST_CASE("FDiscover_Test")
{
    const auto stats = FDiscover::analyze("tests/data/fusion/fusions.out");

    REQUIRE(stats.p.m.sn() == Approx(0.7916666667));
    REQUIRE(stats.p.m.sp() == Approx(1.0));
    REQUIRE(stats.p.m.nr == 24);
}