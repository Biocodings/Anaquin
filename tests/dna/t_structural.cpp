#include <catch.hpp>
#include "dna/structural.hpp"

using namespace Spike;

TEST_CASE("DNA_Variation_Simulation")
{
    const auto r = Structural::analyze("tests/data/dna_sims/DNA.flat.chrT.vcf");
    
    REQUIRE(r.p.m.sp() == Approx(0.0100671141));
    REQUIRE(r.p.m.sn() == Approx(0.012244898));

    REQUIRE(r.p_al.m.sp() == Approx(1.0));
    REQUIRE(r.p_al.m.sn() == Approx(0.9551020408));
    
    REQUIRE(r.p_gt.m.sp() == Approx(1.0));
    REQUIRE(r.p_gt.m.sn() == Approx(0.012244898));

    REQUIRE(r.p_af.m.sp() == Approx(1.0));
    REQUIRE(r.p_af.m.sn() == Approx(0.012244898));
    
    REQUIRE(r.p_l.m.sp() == Approx(1.0));
    REQUIRE(r.p_l.m.sn() == Approx(0.9551020408));
}