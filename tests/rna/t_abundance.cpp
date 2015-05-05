#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_TMap")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/rna.transcripts_an.gtf.tmap");

    REQUIRE(r.lm.r == Approx(0.6833372));
    REQUIRE(r.lm.r == Approx(0.6833372));
    //REQUIRE(r.slope == Approx(0.6833372));
    //REQUIRE(r.r2 == Approx(0.6833372));
}