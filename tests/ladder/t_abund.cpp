#include <catch.hpp>
#include "ladder/l_abund.hpp"

static std::string _output, _error;

extern int parse_options(const std::string &, std::string &, std::string &);

using namespace Anaquin;

TEST_CASE("LAbund_Command")
{
    const int s = parse_options("-t LadderAbund -m data/ladder/CON.v3.mix.csv -usam tests/data/ladder/aligned_A.sam", _output, _error);

    REQUIRE(s == 0);
    REQUIRE(_output.find("Ladder Analysis") != std::string::npos);
}

TEST_CASE("LAbund_Test")
{
    const auto r = LAbund::analyze("tests/data/ladder/aligned_A.sam");

    REQUIRE(r.expTotal == 2785280);
    REQUIRE(r.obsTotal == 834);

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

    REQUIRE(r.expect.at("C_16_A") == Approx(0.0117647059));
    REQUIRE(r.expect.at("C_16_B") == Approx(0.0235294118));
    REQUIRE(r.expect.at("C_16_C") == Approx(0.0470588235));
    REQUIRE(r.expect.at("C_16_D") == Approx(0.0941176471));

    REQUIRE(r.adjusted.at("C_16_A") == Approx(0.0154114804));
    REQUIRE(r.adjusted.at("C_16_B") == Approx(0.0263863226));
    REQUIRE(r.adjusted.at("C_16_C") == Approx(0.0558082398));
    REQUIRE(r.adjusted.at("C_16_D") == Approx(0.0971390282));
}