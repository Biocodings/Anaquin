#include <catch.hpp>
#include "unit/test.hpp"

using namespace Anaquin;

TEST_CASE("Anaquin_Version")
{
    const auto r = Test::test("-v");
    
    REQUIRE(r.status == 0);
    REQUIRE(r.output == "v0.7\n");
}