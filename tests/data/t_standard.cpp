#include <catch.hpp>
#include "unit/test.hpp"
#include "data/standard.hpp"

using namespace Anaquin;

TEST_CASE("isSynthetic_T1")
{
    REQUIRE(Standard::isSynthetic("chrT"));
    REQUIRE(Standard::isSynthetic("chrIS"));
    REQUIRE(!Standard::isSynthetic("chris"));
    REQUIRE(!Standard::isSynthetic("CHRT"));
    REQUIRE(!Standard::isSynthetic("chrt"));
    REQUIRE(!Standard::isSynthetic("Chris"));
    REQUIRE(!Standard::isSynthetic("CHRis"));
    REQUIRE(!Standard::isSynthetic("ChrT"));
}

TEST_CASE("isSynthetic_T2")
{
    REQUIRE(!Standard::isSynthetic("chr1"));
    REQUIRE(!Standard::isSynthetic("chr2"));
    REQUIRE(!Standard::isSynthetic("chr3"));
    REQUIRE(!Standard::isSynthetic("chr4"));
    REQUIRE(!Standard::isSynthetic("chr5"));
    REQUIRE(!Standard::isSynthetic("chr6"));
    REQUIRE(!Standard::isSynthetic("chr7"));
    REQUIRE(!Standard::isSynthetic("chr8"));
}