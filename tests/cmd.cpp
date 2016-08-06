#include <catch.hpp>
#include "unit/test.hpp"

using namespace Anaquin;

TEST_CASE("Anaquin_Version")
{
    const auto r = Test::test("-v");
    
    REQUIRE(r.status == 0);
    REQUIRE(r.output == "v0.99.0\n");
}

TEST_CASE("Anaquin_RnaAlign_HelpShort")
{
    REQUIRE(Test::test("RnaAlign -h").status == 0);
}

TEST_CASE("Anaquin_RnaAlign_HelpLong")
{
    REQUIRE(Test::test("RnaAlign --help").status == 0);
}