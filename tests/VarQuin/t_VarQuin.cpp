#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/VarQuin.hpp"

using namespace Anaquin;

TEST_CASE("isReverseChr")
{
    REQUIRE(isReverseChr("chrev1"));
    REQUIRE(isReverseChr("chrev2"));
    REQUIRE(isReverseChr("chrev3"));
    REQUIRE(isReverseChr("chrevX"));
    
    REQUIRE(isReverseChr("rev1"));
    REQUIRE(isReverseChr("rev2"));
    REQUIRE(isReverseChr("rev3"));
    REQUIRE(isReverseChr("revX"));
    
    REQUIRE(!isReverseChr("chr1"));
    REQUIRE(!isReverseChr("chr2"));
    REQUIRE(!isReverseChr("chr3"));
    REQUIRE(!isReverseChr("chrX"));

    REQUIRE(!isReverseChr("1"));
    REQUIRE(!isReverseChr("2"));
    REQUIRE(!isReverseChr("3"));
    REQUIRE(!isReverseChr("X"));
}
