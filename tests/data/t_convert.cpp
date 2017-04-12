#include <catch.hpp>
#include "tools/tools.hpp"

using namespace Anaquin;

TEST_CASE("ld2ns_Test")
{
    const auto r1 = ld2ss(0.01);
    const auto r2 = ld2ss(0.00000000000000000000000000000001012012012012102);
    const auto r3 = ld2ss(0.99);
    
    REQUIRE(r1 == "1.000000e-02");
    REQUIRE(r2 == "1.012012e-32");
    REQUIRE(r3 == "9.900000e-01");
}

TEST_CASE("ss2ld_Test")
{
    const auto r1 = ss2ld("1e01");
    const auto r2 = ss2ld("1e-01");
    const auto r3 = ss2ld("1e-0300");
    const auto r4 = ss2ld("1e-01");
    
    const auto r5 = ld2ss(r1);
    const auto r6 = ld2ss(r2);
    const auto r7 = ld2ss(r3);
    const auto r8 = ld2ss(r4);
    
    REQUIRE(r1 == 10.0);
    REQUIRE(r2 == Approx(0.1));
    REQUIRE(r3 == Approx(1e-300));

    REQUIRE(r1 > r2);
    REQUIRE(r3 < r2);
    REQUIRE(r4 == r4);
    
    REQUIRE(r5 == "1.000000e+01");
    REQUIRE(r6 == "1.000000e-01");
    REQUIRE(r7 == "1.000000e-300");
    REQUIRE(r8 == "1.000000e-01");
}
