#include <catch.hpp>
#include "aligner.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

TEST_CASE("RNA_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution. It's obviously independent to the standards.
    const auto stats = Aligner::analyze("tests/data/cufflinks_test.sam", Aligner::OutputMode::None);
    
    REQUIRE(isnan(stats.m.sp()));
    REQUIRE(1 == stats.m.sn());
    REQUIRE(0 == stats.m.tp);
    REQUIRE(0 == stats.m.fp);
    REQUIRE(0 == stats.m.fn);
    REQUIRE(3271 == stats.m.tn);
    REQUIRE(0 == stats.nr);
    REQUIRE(3307 == stats.nq);
    REQUIRE(0 == stats.dilution);
}

TEST_CASE("RNA_Simulation_Base")
{
    for (auto ex : exts)
    {
        const auto stats = Aligner::analyze("tests/data/rna_sims/accepted_hits." + ex, Aligner::OutputMode::None);

        REQUIRE(stats.n == 9997);
        REQUIRE(stats.m.tp == 9997);
        REQUIRE(stats.m.fp == 0);
        REQUIRE(stats.m.fn == 0);
        REQUIRE(stats.m.tn == 0);
        REQUIRE(stats.nr == 9997);
        REQUIRE(stats.dilution == 1);
        REQUIRE(stats.nq == 0);
    }
}

TEST_CASE("RNA_Simulation_Splicing")
{
    for (auto ex : exts)
    {
        Aligner::AlignerOptions options;

        options.output = Aligner::OutputMode::None;
        options.mode = Aligner::SpliceAlign;
        
        const auto stats = Aligner::analyze("tests/data/rna_sims/accepted_hits." + ex, options);
        
        REQUIRE(stats.n == 4235);
        REQUIRE(stats.m.sp() == 1);
        REQUIRE(isnan(stats.m.sn()));
        REQUIRE(stats.m.tp == 4235);
        REQUIRE(stats.m.fp == 0);
        REQUIRE(stats.m.fn == 0);
        REQUIRE(stats.m.tn == 0);
        REQUIRE(stats.nr == 4235);
        REQUIRE(stats.dilution == 1);
        REQUIRE(stats.nq == 0);
    }
}

TEST_CASE("RNA_Simulation_Exon")
{
    for (auto ex : exts)
    {
        Aligner::AlignerOptions options;

        options.output = Aligner::OutputMode::None;        
        options.mode = Aligner::ExonAlign;
        
        const auto stats = Aligner::analyze("tests/data/rna_sims/accepted_hits." + ex, options);
        
        REQUIRE(stats.n == 5762);
        REQUIRE(stats.m.tp == 5762);
        REQUIRE(stats.m.fp == 0);
        REQUIRE(stats.m.fn == 0);
        REQUIRE(stats.m.tn == 0);
        REQUIRE(stats.nr == 5762);
        REQUIRE(stats.dilution == 1);
        REQUIRE(stats.nq == 0);
        REQUIRE(stats.m.sp() == 1);
        REQUIRE(isnan(stats.m.sn()));
    }
}