#include <catch.hpp>
#include "MetaQuin/MetaQuin.hpp"

using namespace Anaquin;

TEST_CASE("IsMetaQuin_Test")
{
    REQUIRE(isMetaQuin("MG_03"));
    REQUIRE(isMetaQuin("MG_53"));

    REQUIRE(isMetaQuin("M1_G"));
    REQUIRE(isMetaQuin("M16_G"));
    
    REQUIRE(isMetaQuin("GC_68_1"));
    REQUIRE(isMetaQuin("GC_74_3"));
    
    REQUIRE(isMetaQuin("MH_03"));
    REQUIRE(isMetaQuin("MH_20"));
}
