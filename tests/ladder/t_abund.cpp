#include <catch.hpp>
#include "ladder/l_abund.hpp"

using namespace Anaquin;

TEST_CASE("LAbund_Test")
{
    const auto r = LAbund::analyze("tests/data/ladder/aligned_A.sam");

    REQUIRE(r.expTotal == 4915200);
    REQUIRE(r.actTotal == 834);

    REQUIRE(r.expect.count("C_16_A") == 1);
    REQUIRE(r.expect.count("C_16_B") == 1);
    REQUIRE(r.expect.count("C_16_C") == 1);
    REQUIRE(r.expect.count("C_16_D") == 1);

    REQUIRE(r.measured.count("C_16_A") == 1);
    REQUIRE(r.measured.count("C_16_B") == 1);
    REQUIRE(r.measured.count("C_16_C") == 1);
    REQUIRE(r.measured.count("C_16_D") == 1);

    REQUIRE(r.adjusted.count("C_16_A") == 1);
    REQUIRE(r.adjusted.count("C_16_B") == 1);
    REQUIRE(r.adjusted.count("C_16_C") == 1);
    REQUIRE(r.adjusted.count("C_16_D") == 1);

    REQUIRE(r.measured.at("C_16_A") == 66);
    REQUIRE(r.measured.at("C_16_B") == 113);
    REQUIRE(r.measured.at("C_16_C") == 239);
    REQUIRE(r.measured.at("C_16_D") == 416);

    REQUIRE(r.normalized.at("C_16_A") == Approx(0.0791366906));
    REQUIRE(r.normalized.at("C_16_B") == Approx(0.1354916067));
    REQUIRE(r.normalized.at("C_16_C") == Approx(0.2865707434));
    REQUIRE(r.normalized.at("C_16_D") == Approx(0.4988009592));

    REQUIRE(r.expect.at("C_16_A") == Approx(0.0666666667));
    REQUIRE(r.expect.at("C_16_B") == Approx(0.1333333333));
    REQUIRE(r.expect.at("C_16_C") == Approx(0.2666666667));
    REQUIRE(r.expect.at("C_16_D") == Approx(0.5333333333));

    REQUIRE(r.adjusted.at("C_16_A") == Approx(0.0873317225));
    REQUIRE(r.adjusted.at("C_16_B") == Approx(0.1495224945));
    REQUIRE(r.adjusted.at("C_16_C") == Approx(0.316246692));
    REQUIRE(r.adjusted.at("C_16_D") == Approx(0.5504544932));
}