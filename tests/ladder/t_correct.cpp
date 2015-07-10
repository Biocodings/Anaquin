#include <catch.hpp>
#include "ladder/l_correct.hpp"

using namespace Spike;

TEST_CASE("Ladder_Correct_Test")
{
    const auto r = LCorrect::analyze("tests/data/con/aligned_A.sam");

    REQUIRE(r.expTotal == 450);
    REQUIRE(r.actTotal == 54);

    REQUIRE(r.expect.count("GA116_A") == 1);
    REQUIRE(r.expect.count("GA116_B") == 1);
    REQUIRE(r.expect.count("GA116_C") == 1);
    REQUIRE(r.expect.count("GA116_D") == 1);
    
    REQUIRE(r.abund.count("GA116_A") == 1);
    REQUIRE(r.abund.count("GA116_B") == 1);
    REQUIRE(r.abund.count("GA116_C") == 1);
    REQUIRE(r.abund.count("GA116_D") == 1);

    REQUIRE(r.correct.count("GA116_A") == 1);
    REQUIRE(r.correct.count("GA116_B") == 1);
    REQUIRE(r.correct.count("GA116_C") == 1);
    REQUIRE(r.correct.count("GA116_D") == 1);

    REQUIRE(r.abund.at("GA116_A") == 2);
    REQUIRE(r.abund.at("GA116_B") == 4);
    REQUIRE(r.abund.at("GA116_C") == 7);
    REQUIRE(r.abund.at("GA116_D") == 17);

    REQUIRE(r.actual.at("GA116_A") == Approx(0.037037037));
    REQUIRE(r.actual.at("GA116_B") == Approx(0.0740740741));
    REQUIRE(r.actual.at("GA116_C") == Approx(0.1296296296));
    REQUIRE(r.actual.at("GA116_D") == Approx(0.3148148148));

    REQUIRE(r.expect.at("GA116_A") == Approx(0.0444444444));
    REQUIRE(r.expect.at("GA116_B") == Approx(0.0888888889));
    REQUIRE(r.expect.at("GA116_C") == Approx(0.1777777778));
    REQUIRE(r.expect.at("GA116_D") == Approx(0.3555555556));

    REQUIRE(r.correct.at("GA116_A") == Approx(0.0415537489));
    REQUIRE(r.correct.at("GA116_B") == Approx(0.0831074977));
    REQUIRE(r.correct.at("GA116_C") == Approx(0.145438121));
    REQUIRE(r.correct.at("GA116_D") == Approx(0.3532068654));
}