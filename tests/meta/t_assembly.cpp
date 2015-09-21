#include <catch.hpp>
#include "unit/test.hpp"
#include "meta/m_assembly.hpp"

using namespace Anaquin;

TEST_CASE("DAssembly_M_1005")
{
    Test::meta();
    const auto r = Velvet::analyze<DAsssembly::Stats<Contig>, Contig>("tests/data/M_1005/contigs.fa");
    
    REQUIRE(r.contigs.size() == 453);
    REQUIRE(r.N80 == 549);
    REQUIRE(r.N50 == 649);
    REQUIRE(r.N20 == 1994);
    REQUIRE(r.min == 511);
    REQUIRE(r.max == 1994);
    REQUIRE(r.sum == 7805);
}