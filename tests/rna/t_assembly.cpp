#include <catch.hpp>
#include "assembly.hpp"

using namespace Spike;

TEST_CASE("Assembly_Simulations_Annotated")
{
    const auto r = Assembly::analyze("tests/data/rna_sims/transcripts_an.gtf");

    REQUIRE(r.n  == 1040);
    REQUIRE(r.nr == 1040);
    REQUIRE(r.nq == 0);
}

TEST_CASE("Assembly_Simulations_Denovo")
{
    //const auto r = Assembly::analyze("tests/data/rna_sims/transcripts_dn.gtf");

    //REQUIRE(r.n  == 1040);
    //REQUIRE(r.nr == 1040);
    //REQUIRE(r.nq == 0);
}