#include <catch.hpp>
#include "meta/m_assembly.hpp"

using namespace Spike;

TEST_CASE("DNAssembly_Contigs")
{
    const auto r = Velvet::parse<DNAsssembly::Stats<Contig>, Contig>("tests/data/meta/A/contigs.fa");
    
    REQUIRE(r.contigs.size() == 23299);
    REQUIRE(r.N80 == 596);
    REQUIRE(r.N50 == 634);
    REQUIRE(r.N20 == 777);
    REQUIRE(r.min == 506);
    REQUIRE(r.max == 931);
    REQUIRE(r.sum == 5300);
}