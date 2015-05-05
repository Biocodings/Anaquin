#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_TMap")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/rna.transcripts_an.gtf.tmap");

    REQUIRE(r.lm.r == Approx(0.6833372));
    REQUIRE(r.lm.c == Approx(2427.4274060494458));
    REQUIRE(r.lm.m == Approx(0.024639005789689628));
}