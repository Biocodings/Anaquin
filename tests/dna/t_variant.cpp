#include <catch.hpp>
#include "dna/d_variant.hpp"

using namespace Spike;

TEST_CASE("DNA_Variant_Simulation")
{
    const auto r = DVariant::analyze("tests/data/dna/VARquin.MixA.v1.vcf");

    REQUIRE(r.m.sp()  == Approx(0.977778));
    REQUIRE(r.m.sn()  == Approx(0.347826));
    REQUIRE(r.covered == Approx(0.0553359684));
}