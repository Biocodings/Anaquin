#include <catch.hpp>
#include "Kallisto.hpp"
#include <iostream>

using namespace Anaquin;

TEST_CASE("Kallisto_1")
{
    const auto x = KBuildIndex("tests/data/sequins.fa", 31);
    REQUIRE(!x.empty());
    
    REQUIRE( KQuery___(x, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACC"));
    REQUIRE( KQuery___(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGG"));
    REQUIRE(!KQuery___(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGA"));
    REQUIRE(!KQuery___(x, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACA"));
}

TEST_CASE("Kallisto_2")
{
    const auto f = KHumanFASTA("tests/data/sequins.fa");
    const auto i = KBuildIndex(f, 31);
    REQUIRE(!i.empty());
    
    REQUIRE( KQuery___(i, "ATCTAACGTAAAAACCCTTATTAATTTCATA"));
    REQUIRE( KQuery___(i, "TGAGAGCAGAGCCTGGGTGAGGCAGAATGAAA"));
    REQUIRE(!KQuery___(i, "ATACTTTAATTATTCCCAAAAATGCAATCTA"));
    REQUIRE(!KQuery___(i, "AAAGTAAGACGGAGTGGGTCCGAGACGAGAGT"));
}
