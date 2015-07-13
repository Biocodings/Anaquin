#include <catch.hpp>
#include "meta/m_assembly.hpp"

using namespace Anaquin;

TEST_CASE("DNAssembly_Contigs")
{
    const auto r = Velvet::parse<DNAsssembly::Stats<Contig>, Contig>("tests/data/meta/A/contigs.fa");
    
    REQUIRE(r.contigs.size() == 453);
    REQUIRE(r.N80 == 549);
    REQUIRE(r.N50 == 649);
    REQUIRE(r.N20 == 1994);
    REQUIRE(r.min == 511);
    REQUIRE(r.max == 1994);
    REQUIRE(r.sum == 7805);
}