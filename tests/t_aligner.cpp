#include <catch.hpp>
#include "aligner.hpp"

TEST_CASE("RNA_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution. It's obviously independent to the standards.
    const auto stats = Aligner::analyze("tests/data/cufflinks_test.sam");
    
    /*
     * There shouldn't be any match. Both sensivity and specificity are NAN because the experiement gives no
     * power to detect anything.
     */
    
    REQUIRE(isnan(stats.m.sp()));
    REQUIRE(1 == stats.m.sn());
    REQUIRE(0 == stats.m.tp);
    REQUIRE(0 == stats.m.fp);
    REQUIRE(0 == stats.m.fn);
    REQUIRE(3271 == stats.m.tn);
    REQUIRE(0 == stats.nr);
    REQUIRE(0 == stats.dilution);
    REQUIRE(3307 == stats.nq);
}

TEST_CASE("RNA_Simulation_Base")
{
    /*
     * Since the SAM file comes from a simulation, it's not surprising that sensitivity is perfect.
     */
    
    const auto stats = Aligner::analyze("tests/data/rna_sims/accepted_hits.sam");

    REQUIRE(stats.n == 9997);
    REQUIRE(stats.m.sp() == 1);
    REQUIRE(isnan(stats.m.sn()));
    REQUIRE(stats.m.tp == 9997);
    REQUIRE(stats.m.fp == 0);
    REQUIRE(stats.m.fn == 0);
    REQUIRE(stats.m.tn == 0);
    REQUIRE(stats.nr == 9997);
    REQUIRE(stats.dilution == 1);
    REQUIRE(stats.nq == 0);
}

TEST_CASE("RNA_Simulation_Splicing")
{
    /*
     * Since the SAM file comes from a simulation of the reference chromosome, it's not surprising that
     * the sensitivity is 100% (from the first 1000 reads).
     */
    
    //const auto stats = Aligner::spliced("tests/data/rna_sims/accepted_hits.sam", AlignerOptions(1000));

    //REQUIRE(1 == stats.m.sp());
    //REQUIRE(0 == stats.m.sn());
}