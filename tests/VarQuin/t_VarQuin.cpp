#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/VarQuin.hpp"

using namespace Anaquin;

TEST_CASE("isRevChr")
{
    REQUIRE(isRevChr("chrev1"));
    REQUIRE(isRevChr("chrev2"));
    REQUIRE(isRevChr("chrev3"));
    REQUIRE(isRevChr("chrevX"));
    
    REQUIRE(isRevChr("rev1"));
    REQUIRE(isRevChr("rev2"));
    REQUIRE(isRevChr("rev3"));
    REQUIRE(isRevChr("revX"));
    
    REQUIRE(!isRevChr("chr1"));
    REQUIRE(!isRevChr("chr2"));
    REQUIRE(!isRevChr("chr3"));
    REQUIRE(!isRevChr("chrX"));

    REQUIRE(!isRevChr("1"));
    REQUIRE(!isRevChr("2"));
    REQUIRE(!isRevChr("3"));
    REQUIRE(!isRevChr("X"));
}
