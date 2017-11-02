#include <catch.hpp>
#include "Kallisto.hpp"
#include <iostream>

using namespace Anaquin;

TEST_CASE("Kallisto_1")
{
    const auto x = KBuildIndex("tests/data/A.V.28.fa", 31);
    REQUIRE(!x.empty());
    
    REQUIRE( KQuery___(x, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACC"));
    REQUIRE( KQuery___(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGG"));
    REQUIRE(!KQuery___(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGA"));
    REQUIRE(!KQuery___(x, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACA"));
}

TEST_CASE("Kallisto_2")
{
    const auto f = KHumanFASTA("tests/data/A.V.28.fa");
    const auto x = KBuildIndex(f, 31);
    
    REQUIRE( KQuery___(x, "GGAAGGGGCAGATTTCGGGGTCTAGGCTTGG"));
    REQUIRE(!KQuery___(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGG"));
}
