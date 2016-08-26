#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_flip.hpp"

using namespace Anaquin;

TEST_CASE("VFlip_IsReverse_1")
{
    const auto x = std::set<ReadName> { "GS_016_V-GV_hom_vars.sim_reads7478" };
    
    REQUIRE(VFlip::isReverse(x, "GS_016_V-GV_hom_vars.sim_reads7478"));
    REQUIRE(VFlip::isReverse(x, "GS_016_V-GV_hom_vars.sim_reads7478/1"));
    REQUIRE(VFlip::isReverse(x, "GS_016_V-GV_hom_vars.sim_reads7478/2"));
    REQUIRE(VFlip::isReverse(x, "@GS_016_V-GV_hom_vars.sim_reads7478/1"));
    REQUIRE(VFlip::isReverse(x, "@GS_016_V-GV_hom_vars.sim_reads7478/2"));
}