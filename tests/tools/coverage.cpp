#include <catch.hpp>
#include "unit/test.hpp"

using namespace Anaquin;

TEST_CASE("Test_Version")
{
    const auto r = Test::test("-v");
    
    REQUIRE(r.status == 0);
    REQUIRE(r.output == "Anaquin v1.1.01\n");
}