#include <catch.hpp>
#include "meta/m_assembly.hpp"

using namespace Spike;

/*
 * Reference: http://quast.bioinf.spbau.ru/manual.html
 */

TEST_CASE("MAssembly_Contigs_2")
{
    const auto r = MAssembly::analyze("tests/data/meta/contigs_2.fa");

    REQUIRE(r.N50  == 594);
    REQUIRE(r.N80  == 594);
    REQUIRE(r.N20  == 594);
    REQUIRE(r.mean == 594);
    REQUIRE(r.min  == 594);
    REQUIRE(r.max  == 594);
    REQUIRE(r.sum  == 594);
    REQUIRE(r.total == 168340);
    REQUIRE(r.contigs.size() == 2764);
}

TEST_CASE("DNAssembly_Contigs_1")
{
    const auto r = Velvet::parse<DNAsssembly::Stats<Contig>, Contig>("tests/data/meta/contigs_1.fa");
    
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