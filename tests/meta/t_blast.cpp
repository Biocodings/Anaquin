#include <catch.hpp>
#include "meta/m_blast.hpp"

using namespace Spike;

TEST_CASE("MBlast_Empty")
{
    REQUIRE_THROWS(MBlast::analyze("tests/data/empty.psl"));
}