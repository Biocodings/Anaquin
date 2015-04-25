#include <catch.hpp>
#include "r_assembly.hpp"

using namespace Spike;

TEST_CASE("Assembly_Simulations_Annotated")
{
    const auto r = RAssembly::analyze("tests/data/rna_sims/transcripts_an.gtf");

    REQUIRE(r.me.nq() == 490);
    
    REQUIRE(r.se.id == "R_1_4_R");
    REQUIRE(r.se.counts == 1);
    REQUIRE(r.se.abund == 152);
    
    REQUIRE(r.st.id == "R_9_2_R");
    REQUIRE(r.st.counts == 1);
    REQUIRE(r.st.abund == 9.25);
    
    
//    REQUIRE(r.n()  == 1040);
  //  REQUIRE(r.nr == 1040);
    //REQUIRE(r.nq == 0);
}

TEST_CASE("Assembly_Simulations_Denovo")
{
    //const auto r = RAssembly::analyze("tests/data/rna_sims/transcripts_dn.gtf");

    //REQUIRE(r.n  == 1040);
    //REQUIRE(r.nr == 1040);
    //REQUIRE(r.nq == 0);
}