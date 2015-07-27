#include <catch.hpp>
#include "fusion/f_express.hpp"

static std::string _output, _error;

extern int parse_options(const std::string &, std::string &, std::string &);

using namespace Anaquin;

TEST_CASE("FExpress_Test")
{
    const auto stats = FExpress::analyze("tests/data/fusion/fusions.out");

    // The linear model associated with the expression
    const auto lm = stats.linear();

    REQUIRE(lm.r  == Approx(0.9827924165));
    REQUIRE(lm.r2 == Approx(0.96387393));
    REQUIRE(lm.m  == Approx(0.8617847724));
    REQUIRE(lm.c  == Approx(3.2714500279));
}

TEST_CASE("FExpress_Command")
{
    const auto cmd = "-t FusionExpress -rfus data/fusion/FUS.v1.ref -m data/fusion/FUS.v3.csv -uout tests/data/fusion/fusions.out";

    const auto status = parse_options(cmd, _output, _error);
    
    REQUIRE(status == 0);
}