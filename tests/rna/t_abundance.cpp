#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_Genes_Tracking")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/genes.fpkm_tracking");

    REQUIRE(r.lm.r2 == Approx(0.9610729944));
    REQUIRE(r.lm.r  == Approx(0.9811371334));
    REQUIRE(r.lm.c  == Approx(-1.9163438325));
    REQUIRE(r.lm.m  == Approx(0.9556247487));

    REQUIRE(r.s.id == "R_9_2");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == 37.0);
}

TEST_CASE("RAbundance_Simulations_Isoforms_Tracking")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/isoforms.fpkm_tracking");
    
    REQUIRE(r.lm.r  == Approx(0.8079915259));
    REQUIRE(r.lm.c  == Approx(-2.7791382851));
    REQUIRE(r.lm.m  == Approx(0.9945833137));
    REQUIRE(r.lm.r2 == Approx(0.6456180206));

    REQUIRE(r.s.id == "R_9_2_R");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == 9.25);
}

TEST_CASE("RAbundance_Simulations_TMap")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/rna.transcripts.gtf.tmap");

    REQUIRE(r.lm.r  == Approx(0.8062709437));
    REQUIRE(r.lm.c  == Approx(-2.7899650615));
    REQUIRE(r.lm.m  == Approx(0.995560013));
    REQUIRE(r.lm.r2 == Approx(0.6427826854));

    REQUIRE(r.s.id == "R_9_2_R");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == 9.25);
}