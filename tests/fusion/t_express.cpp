#include <catch.hpp>
#include "fusion/f_express.hpp"

using namespace Anaquin;

TEST_CASE("FExpress_Simulated")
{
    //./anaquin -t FusionExpress -uout genes.fpkm_tracking
    FExpress::analyze("/Users/tedwong/Sources/QA/genes.fpkm_tracking");
}