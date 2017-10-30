#include <catch.hpp>
#include "Kallisto.hpp"

using namespace Anaquin;

TEST_CASE("Kallisto_Index")
{
    const auto index = KBuildIndex("tests/data/sequins.fa", 31);
    REQUIRE(!index.empty());
    
    REQUIRE(KQuery(index, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACC"));
    REQUIRE(KQuery(index, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGG"));
    REQUIRE(!KQuery(index, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGA"));
    REQUIRE(!KQuery(index, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACA"));
}
