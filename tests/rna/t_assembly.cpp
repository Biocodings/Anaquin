#include <catch.hpp>
#include "r_assembly.hpp"

using namespace Spike;

TEST_CASE("Assembly_Simulations_All_Filtered")
{
    RAssembly::Options opts;
    const auto &s = Standard::instance();

    for (auto i: s.r_seqs_iA)
    {
        opts.filters.insert(i.first);
    }
    
    const auto r = RAssembly::analyze("tests/data/rna_sims/transcripts_an.gtf", opts);
    
    REQUIRE(r.me.nq() == 0);
    REQUIRE(r.me.nr() == 498);
    REQUIRE(r.mt.nq() == 0);
    REQUIRE(r.mt.nr() == 61);

    REQUIRE(r.me.tp() == 0);
    REQUIRE(r.me.fn() == 498);
}

TEST_CASE("Assembly_Simulations_Annotated")
{
    const auto r = RAssembly::analyze("tests/data/rna_sims/transcripts_an.gtf");

    REQUIRE(r.me.nq() == 490);
    REQUIRE(r.me.nr() == 498);
    REQUIRE(r.mt.nq() == 61);
    REQUIRE(r.mt.nr() == 61);
    
    REQUIRE(r.me.tp() == 490);
    REQUIRE(r.me.fn() == 8);
    
    REQUIRE(r.se.id == "R_1_4_R");
    REQUIRE(r.se.counts == 1);
    REQUIRE(r.se.abund == 152);
    
    REQUIRE(r.st.id == "R_9_2_R");
    REQUIRE(r.st.counts == 1);
    REQUIRE(r.st.abund == 9.25);
}

TEST_CASE("Assembly_Simulations_Denovo")
{
    int a = 1;
    a = 1;
    //const auto r = RAssembly::analyze("tests/data/rna_sims/transcripts_dn.gtf");

    //REQUIRE(r.n  == 1040);
    //REQUIRE(r.nr == 1040);
    //REQUIRE(r.nq == 0);
}