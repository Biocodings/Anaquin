#include <catch.hpp>
#include "test.hpp"

using namespace Anaquin;

TEST_CASE("Anaquin_RnaAlign_HelpShort")
{
    REQUIRE(runTest("RnaAlign -h").status == 0);
}

TEST_CASE("Anaquin_RnaAlign_HelpLong")
{
    REQUIRE(runTest("RnaAlign --help").status == 0);
}
