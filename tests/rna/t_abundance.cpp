#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_Genes_Tracking")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/genes.fpkm_tracking");
    
    REQUIRE(r.lm.r == Approx(0.7960605027));
    REQUIRE(r.lm.c == Approx(104.1467805074));
    REQUIRE(r.lm.m == Approx(0.0264169601));
    
    REQUIRE(r.s.id == "R_9_2");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == 37.0);
}

TEST_CASE("RAbundance_Simulations_Isoforms_Tracking")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/isoforms.fpkm_tracking");
    
    REQUIRE(r.lm.r == Approx(0.6833372103));
    REQUIRE(r.lm.c == Approx(2427.4335867495));
    REQUIRE(r.lm.m == Approx(0.0246389856));

    REQUIRE(r.s.id == "R_9_2_R");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == 9.25);
}

TEST_CASE("RAbundance_Simulations_TMap")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/rna.transcripts_an.gtf.tmap");

    REQUIRE(r.lm.r == Approx(0.6833372));
    REQUIRE(r.lm.c == Approx(2427.4274060494458));
    REQUIRE(r.lm.m == Approx(0.024639005789689628));

    REQUIRE(r.s.id == "R_9_2_R");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == 9.25);
}