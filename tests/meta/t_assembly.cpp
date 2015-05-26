#include <catch.hpp>
#include "meta/m_assembly.hpp"

using namespace Spike;

TEST_CASE("DNAssembly_Contigs")
{
    const auto r = DNAsssembly::stats("tests/data/meta/contigs.fa");

    REQUIRE(r.contigs.size() == 63);

    REQUIRE(r.contigs[0].id == "NODE_1_length_4075_cov_20.748220");
    REQUIRE(r.contigs[1].id == "NODE_2_length_1635_cov_21.235474");
    REQUIRE(r.contigs[2].id == "NODE_3_length_2338_cov_20.628742");
    REQUIRE(r.contigs[3].id == "NODE_4_length_1996_cov_19.849699");

    REQUIRE(r.N80  == 1836);
    REQUIRE(r.N50  == 2846);
    REQUIRE(r.N20  == 4703);
    REQUIRE(r.mean == 3585);
    REQUIRE(r.min  == 67);
    REQUIRE(r.max  == 9109);
    REQUIRE(r.sum  == 139916);
}

TEST_CASE("MAssembly_Contigs")
{
    const auto r = MAssembly::analyze("tests/data/meta/contigs.fa");

    REQUIRE(r.dstats.contigs.size() == 63);
    
    REQUIRE(r.dstats.contigs[0].id == "NODE_1_length_4075_cov_20.748220");
    REQUIRE(r.dstats.contigs[1].id == "NODE_2_length_1635_cov_21.235474");
    REQUIRE(r.dstats.contigs[2].id == "NODE_3_length_2338_cov_20.628742");
    REQUIRE(r.dstats.contigs[3].id == "NODE_4_length_1996_cov_19.849699");
    
    REQUIRE(r.dstats.N80  == 1836);
    REQUIRE(r.dstats.N50  == 2846);
    REQUIRE(r.dstats.N20  == 4703);
    REQUIRE(r.dstats.mean == 3585);
    REQUIRE(r.dstats.min  == 67);
    REQUIRE(r.dstats.max  == 9109);
    REQUIRE(r.dstats.sum  == 139916);
}