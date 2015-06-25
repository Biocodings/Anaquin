#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_Isoforms_A1")
{
    const auto r  = RAbundance::analyze("tests/data/rna/A1/isoforms.fpkm_tracking");
    const auto lm = r.linear();

    REQUIRE(lm.r  == Approx(0.8143944753));
    REQUIRE(lm.m  == Approx(0.9300431127));
    REQUIRE(lm.r2 == Approx(0.6591315122));

    REQUIRE(r.s.id == "R2_38_1");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(1.0));
}

TEST_CASE("RAbundance_Simulations_Genes_A1")
{
    const auto r  = RAbundance::analyze("tests/data/rna/A1/genes.fpkm_tracking");
    const auto lm = r.linear();
    
    REQUIRE(lm.r  == Approx(0.8886262976));
    REQUIRE(lm.m  == Approx(0.9453642698));
    REQUIRE(lm.r2 == Approx(0.7841213466));
    
    REQUIRE(r.s.id == "R2_33");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(7.0));
}