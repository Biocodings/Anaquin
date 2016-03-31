#include <catch.hpp>
#include "unit/test.hpp"
#include "MetaQuin/m_blat.hpp"

using namespace Anaquin;

TEST_CASE("MBlast_Empty")
{
    Test::meta();
    REQUIRE_THROWS(MBlat::report("tests/examples/empty.psl"));
}
