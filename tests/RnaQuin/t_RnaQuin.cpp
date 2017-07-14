#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/RnaQuin.hpp"

using namespace Anaquin;

TEST_CASE("isRNARevChr_1")
{
    REQUIRE(isRNARevChr("chrIS"));
    REQUIRE(!isRNARevChr("chrT"));
    REQUIRE(!isRNARevChr("chris"));
    REQUIRE(!isRNARevChr("CHRT"));
    REQUIRE(!isRNARevChr("chrt"));
    REQUIRE(!isRNARevChr("Chris"));
    REQUIRE(!isRNARevChr("CHRis"));
    REQUIRE(!isRNARevChr("ChrT"));
}

TEST_CASE("isRNARevChr_2")
{
    REQUIRE(!isRNARevChr("chr1"));
    REQUIRE(!isRNARevChr("chr2"));
    REQUIRE(!isRNARevChr("chr3"));
    REQUIRE(!isRNARevChr("chr4"));
    REQUIRE(!isRNARevChr("chr5"));
    REQUIRE(!isRNARevChr("chr6"));
    REQUIRE(!isRNARevChr("chr7"));
    REQUIRE(!isRNARevChr("chr8"));
}
