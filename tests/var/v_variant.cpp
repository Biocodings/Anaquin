#include <catch.hpp>
#include "var/v_variant.hpp"

using namespace Anaquin;

TEST_CASE("Variant_Variant_Simulation")
{
    const auto r  = DVariant::analyze("tests/data/dna/VARquin.MixA.v1.vcf");
    const auto lm = r.linear();

    REQUIRE(r.m.sp()  == Approx(0.977778));
    REQUIRE(r.m.sn()  == Approx(0.347826));
    REQUIRE(r.covered == Approx(0.0553359684));

    REQUIRE(lm.r  == Approx(0.954813));
    REQUIRE(lm.r2 == Approx(0.910642));
    REQUIRE(lm.m  == Approx(0.925754));
}