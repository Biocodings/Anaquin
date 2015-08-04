#include <catch.hpp>
#include "fusion/f_express.hpp"

using namespace Anaquin;

TEST_CASE("FExpress_Simulated")
{
    const auto r = FExpress::analyze("tests/data/fusion/simulated/normal/genes.fpkm_tracking");
    
    // The linear model associated with the expression
    const auto lm = r.linear();

    
    
    //./anaquin -t FusionExpress -ugtrack genes.fpkm_tracking
    //FExpress::analyze("/Users/tedwong/Sources/QA/genes.fpkm_tracking");
}