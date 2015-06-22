#include <catch.hpp>
#include "meta/m_assembly.hpp"

using namespace Spike;

TEST_CASE("DNAssembly_Contigs")
{
    const auto r = Velvet::parse<DNAsssembly::Stats<Contig>, Contig>("tests/data/meta/A/contigs.fa");
    
    REQUIRE(r.contigs.size() == 63);

    REQUIRE(r.contigs.count("NODE_1_length_4075_cov_20.748220"));
    REQUIRE(r.contigs.count("NODE_2_length_1635_cov_21.235474"));
    REQUIRE(r.contigs.count("NODE_3_length_2338_cov_20.628742"));
    REQUIRE(r.contigs.count("NODE_4_length_1996_cov_19.849699"));

    REQUIRE(r.N80  == 1836);
    REQUIRE(r.N50  == 2846);
    REQUIRE(r.N20  == 4703);
    REQUIRE(r.mean == 3610);
    REQUIRE(r.min  == 508);
    REQUIRE(r.max  == 9109);
    REQUIRE(r.sum  == 138888);
}