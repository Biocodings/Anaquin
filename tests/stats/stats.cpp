#include <catch.hpp>
#include "stats/stats.hpp"

using namespace Anaquin;

TEST_CASE("SReals_Summary")
{
    SReals x;
    
    x.add(1.00);
    x.add(2.00);
    x.add(3.00);
    
    REQUIRE(STRING(x) == "2.00 ± 1.00");
}

TEST_CASE("SStrings_Summary")
{
    SStrings x;
    
    x.add("1");
    x.add("2");
    x.add("3");
    
    REQUIRE(STRING(x) == "1, 2, 3");
}