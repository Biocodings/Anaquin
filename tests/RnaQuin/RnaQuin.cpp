#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/RnaQuin.hpp"

using namespace Anaquin;

TEST_CASE("isRnaQuin_1")
{
    REQUIRE(RnaQuin::isRnaQuin("chrIS"));
    REQUIRE(!RnaQuin::isRnaQuin("chrT"));
    REQUIRE(!RnaQuin::isRnaQuin("chris"));
    REQUIRE(!RnaQuin::isRnaQuin("CHRT"));
    REQUIRE(!RnaQuin::isRnaQuin("chrt"));
    REQUIRE(!RnaQuin::isRnaQuin("Chris"));
    REQUIRE(!RnaQuin::isRnaQuin("CHRis"));
    REQUIRE(!RnaQuin::isRnaQuin("ChrT"));
}

TEST_CASE("isRnaQuin_2")
{
    REQUIRE(!RnaQuin::isRnaQuin("chr1"));
    REQUIRE(!RnaQuin::isRnaQuin("chr2"));
    REQUIRE(!RnaQuin::isRnaQuin("chr3"));
    REQUIRE(!RnaQuin::isRnaQuin("chr4"));
    REQUIRE(!RnaQuin::isRnaQuin("chr5"));
    REQUIRE(!RnaQuin::isRnaQuin("chr6"));
    REQUIRE(!RnaQuin::isRnaQuin("chr7"));
    REQUIRE(!RnaQuin::isRnaQuin("chr8"));
}
