#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_Isoforms_A1")
{
    const auto r  = RAbundance::analyze("tests/data/rna/A1/isoforms.fpkm_tracking");
    const auto lm = r.linear();

    REQUIRE(lm.r  == Approx(0.8315693535));
    REQUIRE(lm.m  == Approx(0.8079902611));
    REQUIRE(lm.r2 == Approx(0.6855750429));

    REQUIRE(r.s.id == "R2_38_1");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(1.0));
}

TEST_CASE("RAbundance_Simulations_Genes_A1")
{
    const auto r  = RAbundance::analyze("tests/data/rna/A1/genes.fpkm_tracking");
    const auto lm = r.linear();
    
    REQUIRE(lm.r  == Approx(0.8934756622));
    REQUIRE(lm.m  == Approx(0.8886562013));
    REQUIRE(lm.r2 == Approx(0.7848520096));
    
    REQUIRE(r.s.id == "R2_33");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(7.0));
}