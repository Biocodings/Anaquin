#include <catch.hpp>
#include "meta/m_blast.hpp"

using namespace Anaquin;

TEST_CASE("MBlast_Empty")
{
    REQUIRE_THROWS(MBlast::analyze("tests/data/empty.psl"));
}