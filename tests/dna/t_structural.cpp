#include <catch.hpp>
#include "dna/structural.hpp"

using namespace Spike;

TEST_CASE("DNA_Variation_Simulation")
{
    const auto r = Structural::analyze("tests/data/dna_sims/DNA.flat.chrT.vcf");
    
    REQUIRE(r.m.sp() == Approx(0.7852348993));
    REQUIRE(r.m.sn() == Approx(0.9551020408));

    REQUIRE(r.ml.sp() == Approx(1.0));
    REQUIRE(r.ml.sn() == Approx(0.9551020408));
}