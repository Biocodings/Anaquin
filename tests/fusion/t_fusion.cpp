#include <catch.hpp>
#include "fusion/f_align.hpp"

using namespace Anaquin;

TEST_CASE("FAlign_Test")
{
    const auto stats = FAlign::analyze("tests/data/fusion/fusions.out");

    REQUIRE(stats.p.m.sn() == Approx(0.7916666667));
    REQUIRE(stats.p.m.sp() == Approx(1.0));
    REQUIRE(stats.p.m.nr == 24);
}