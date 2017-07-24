#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/RnaQuin.hpp"

using namespace Anaquin;

TEST_CASE("isChrIS_1")
{
    REQUIRE(isChrIS("chrIS"));
    REQUIRE(!isChrIS("chrT"));
    REQUIRE(!isChrIS("chris"));
    REQUIRE(!isChrIS("CHRT"));
    REQUIRE(!isChrIS("chrt"));
    REQUIRE(!isChrIS("Chris"));
    REQUIRE(!isChrIS("CHRis"));
    REQUIRE(!isChrIS("ChrT"));
}

TEST_CASE("isChrIS_2")
{
    REQUIRE(!isChrIS("chr1"));
    REQUIRE(!isChrIS("chr2"));
    REQUIRE(!isChrIS("chr3"));
    REQUIRE(!isChrIS("chr4"));
    REQUIRE(!isChrIS("chr5"));
    REQUIRE(!isChrIS("chr6"));
    REQUIRE(!isChrIS("chr7"));
    REQUIRE(!isChrIS("chr8"));
}
