#include <catch.hpp>
#include "fusion/f_express.hpp"

using namespace Anaquin;

TEST_CASE("FExpress_Simulated")
{
    const auto r = FExpress::analyze("tests/data/fusion/simulated/normal/genes.fpkm_tracking");
    
    // The linear model associated with the expression
    const auto lm = r.linear();

    REQUIRE(lm.m  == Approx(0.67651639834586408));
    REQUIRE(lm.r2 == Approx(0.085666237321215699));
    REQUIRE(lm.r  == Approx(0.35414669160789208));
}