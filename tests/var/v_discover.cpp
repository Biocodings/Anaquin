#include <catch.hpp>
#include "var/v_discover.hpp"

using namespace Anaquin;

TEST_CASE("VarDiscover_Test")
{
    const auto r  = VDiscover::analyze("tests/data/var/VARquin.MixA.v1.vcf");
    const auto lm = r.linear();

    REQUIRE(r.m.sp()  == Approx(0.977778));
    REQUIRE(r.m.sn()  == Approx(0.347826));
    //REQUIRE(r.covered == Approx(0.0553359684));
}