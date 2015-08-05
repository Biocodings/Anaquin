#include <catch.hpp>
#include "unit/test.hpp"
#include "meta/m_blast.hpp"

using namespace Anaquin;

TEST_CASE("MBlast_Empty")
{
    Test::meta();
    REQUIRE_THROWS(MBlast::analyze("tests/data/empty.psl"));
}