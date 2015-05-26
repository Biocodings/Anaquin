#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("RAbundance_Simulations_Isoforms_Tracking")
{
    const auto r = RAbundance::analyze("/Users/tedwong/Sources/QA/isoforms.fpkm_tracking");
    
    REQUIRE(r.lm.r  == Approx(0.7835626667));
    REQUIRE(r.lm.c  == Approx(3.5085392148));
    REQUIRE(r.lm.m  == Approx(0.9596871734));
    REQUIRE(r.lm.r2 == Approx(0.6059281704));

    REQUIRE(r.s.id == "R_9_2_R");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(0.0190734863));
}

TEST_CASE("RAbundance_Simulations_TMap")
{
    const auto r = RAbundance::analyze("tests/data/rna/rna.transcripts.gtf.tmap");

    REQUIRE(r.lm.r  == Approx(0.7817439443));
    REQUIRE(r.lm.c  == Approx(3.5041984197));
    REQUIRE(r.lm.m  == Approx(0.960445118));
    REQUIRE(r.lm.r2 == Approx(0.6030220026));

    REQUIRE(r.s.id == "R_9_2_R");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(0.0190734863));
}