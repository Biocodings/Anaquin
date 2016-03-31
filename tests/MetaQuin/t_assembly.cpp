#include <catch.hpp>
#include "unit/test.hpp"
#include "MetaQuin/m_assembly.hpp"

using namespace Anaquin;

TEST_CASE("DAssembly_M_1005")
{
    Test::meta();
    const auto r = Velvet::analyze<DAsssembly::Stats<Contig>, Contig>("tests/data/M_1005/contigs.fa");
    
    REQUIRE(r.contigs.size() == 453);
    REQUIRE(r.N80 == 61);
    REQUIRE(r.N50 == 153);
    REQUIRE(r.N20 == 440);
    REQUIRE(r.min == 61);
    REQUIRE(r.max == 1994);
    REQUIRE(r.sum == 57338);
}
