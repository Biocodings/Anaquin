#include <catch.hpp>
#include "rna/aligner.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

TEST_CASE("Align_RNA_Sims_2")
{
    for (auto ex : exts)
    {
        const auto r = Aligner::analyze("tests/data/rna_sims_2/align/accepted_hits." + ex);

        REQUIRE(r.n == 98041);
        REQUIRE(r.sens.id == "R_5_2");
        REQUIRE(r.sens.counts == 15);
        REQUIRE(r.sens.exp == 9765.0);
    }
}

TEST_CASE("Align_RNA_Sims_Exon")
{
    for (auto ex : exts)
    {
        Aligner::Options options;
        options.mode = Aligner::AlignExon;
        const auto r = Aligner::analyze("tests/data/rna_sims/align/accepted_hits." + ex, options);
        
        REQUIRE(r.n == 5762);
        REQUIRE(r.m.tp == 5762);
        REQUIRE(r.m.fp == 0);
        REQUIRE(r.m.fn == 0);
        REQUIRE(r.m.tn == 0);
        REQUIRE(r.nr == 5762);
        REQUIRE(r.dilution() == 1);
        REQUIRE(r.nq == 0);
        REQUIRE(r.m.sp() == 1);
        REQUIRE(isnan(r.m.sn()));
    }
}

TEST_CASE("Align_RNA_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution. It's obviously independent.
    const auto r = Aligner::analyze("tests/data/cufflinks_test.sam");

    REQUIRE(isnan(r.m.sp()));
    REQUIRE(1 == r.m.sn());
    REQUIRE(0 == r.m.tp);
    REQUIRE(0 == r.m.fp);
    REQUIRE(0 == r.m.fn);
    REQUIRE(3271 == r.m.tn);
    REQUIRE(0 == r.nr);
    REQUIRE(3307 == r.nq);
    REQUIRE(0 == r.dilution());
}

TEST_CASE("Align_RNA_Sims_Base")
{
    for (auto ex : exts)
    {
        const auto r = Aligner::analyze("tests/data/rna_sims/align/accepted_hits." + ex);

        REQUIRE(r.n == 9997);
        REQUIRE(r.m.tp == 9997);
        REQUIRE(r.m.fp == 0);
        REQUIRE(r.m.fn == 0);
        REQUIRE(r.m.tn == 0);
        REQUIRE(r.nr == 9997);
        REQUIRE(r.dilution() == 1);
        REQUIRE(r.nq == 0);
    }
}

TEST_CASE("Align_RNA_Sims_Splicing")
{
    for (auto ex : exts)
    {
        Aligner::Options options;
        options.mode = Aligner::AlignSplice;
        const auto r = Aligner::analyze("tests/data/rna_sims/align/accepted_hits." + ex, options);
        
        REQUIRE(r.n == 4235);
        REQUIRE(r.m.sp() == 1);
        REQUIRE(isnan(r.m.sn()));
        REQUIRE(r.m.tp == 4235);
        REQUIRE(r.m.fp == 0);
        REQUIRE(r.m.fn == 0);
        REQUIRE(r.m.tn == 0);
        REQUIRE(r.nr == 4235);
        REQUIRE(r.dilution() == 1);
        REQUIRE(r.nq == 0);
    }
}