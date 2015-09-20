#include <catch.hpp>
#include "unit/test.hpp"
#include "fusion/f_express.hpp"

using namespace Anaquin;

TEST_CASE("FExpress_F_1000")
{
    Test::fusionA();
    
    const auto r = Test::test("-t FusionExpress -m data/fusion/MFU007.v013.csv -rfus data/fusion/AFU004.v032.ref -soft star -uout tests/data/F_1000/star-fusion.fusion_candidates.txt");
    
    REQUIRE(r.status == 0);
    
    Test::fusionA();
    
    const auto stats = FExpress::report("tests/data/F_1000/star-fusion.fusion_candidates.txt",
                                         FDiscover::Options(FDiscover::Software::Star));
    const auto lm = stats.linear();
    
    REQUIRE(lm.r  == Approx(0.9729157505));
    REQUIRE(lm.m  == Approx(0.9770652147));
    REQUIRE(lm.r2 == Approx(0.9465650576));
}