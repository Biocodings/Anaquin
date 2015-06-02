#include <catch.hpp>
#include "meta/m_assembly.hpp"

using namespace Spike;

TEST_CASE("dsdsaklskladlksad")
{
    const auto r = MAssembly::analyze("/Users/tedwong/Sources/QA/contigs.fa");
    
    
    
    REQUIRE(r.ds.N50  == 594);
    REQUIRE(r.ds.N80  == 594);
    REQUIRE(r.ds.N20  == 594);
    REQUIRE(r.ds.mean == 594);
    REQUIRE(r.ds.min  == 594);
    REQUIRE(r.ds.max  == 594);
    REQUIRE(r.ds.sum  == 594);
}


TEST_CASE("MAssembly_Contigs_2")
{
    const auto r = MAssembly::analyze("tests/data/meta/contigs_2.fa");
    
    REQUIRE(r.ds.N50  == 594);
    REQUIRE(r.ds.N80  == 594);
    REQUIRE(r.ds.N20  == 594);
    REQUIRE(r.ds.mean == 594);
    REQUIRE(r.ds.min  == 594);
    REQUIRE(r.ds.max  == 594);
    REQUIRE(r.ds.sum  == 594);
}

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
    REQUIRE(r.mean == 3610);
    REQUIRE(r.min  == 508);
    REQUIRE(r.max  == 9109);
    REQUIRE(r.sum  == 138888);
}
