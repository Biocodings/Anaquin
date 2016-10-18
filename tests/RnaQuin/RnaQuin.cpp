#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/RnaQuin.hpp"

using namespace Anaquin;

TEST_CASE("isRnaQuin_1")
{
    REQUIRE(isRnaQuin("chrIS"));
    REQUIRE(!isRnaQuin("chrT"));
    REQUIRE(!isRnaQuin("chris"));
    REQUIRE(!isRnaQuin("CHRT"));
    REQUIRE(!isRnaQuin("chrt"));
    REQUIRE(!isRnaQuin("Chris"));
    REQUIRE(!isRnaQuin("CHRis"));
    REQUIRE(!isRnaQuin("ChrT"));
}

TEST_CASE("isRnaQuin_2")
{
    REQUIRE(!isRnaQuin("chr1"));
    REQUIRE(!isRnaQuin("chr2"));
    REQUIRE(!isRnaQuin("chr3"));
    REQUIRE(!isRnaQuin("chr4"));
    REQUIRE(!isRnaQuin("chr5"));
    REQUIRE(!isRnaQuin("chr6"));
    REQUIRE(!isRnaQuin("chr7"));
    REQUIRE(!isRnaQuin("chr8"));
}
