#include <catch.hpp>
#include "fusion/f_express.hpp"

static std::string _output, _error;

extern int parse_options(const std::string &, std::string &, std::string &);

using namespace Anaquin;

TEST_CASE("FExpress_10K")
{
    const auto cmd = "-t FusionExpress -rfus tests/data/fusion/10K/FUS.v1.ref -m tests/data/fusion/10K/FUSE_mixtures_v3.csv -uout tests/data/fusion/10K/fusions.out";
    
    const auto status = parse_options(cmd, _output, _error);

    REQUIRE(status == 0);
}

TEST_CASE("FExpress_Simulated")
{
    const auto stats = FExpress::analyze("tests/data/fusion/fusions.out");

    // The linear model associated with the expression
    const auto lm = stats.linear();

    REQUIRE(lm.r  == Approx(0.9489461887));
    REQUIRE(lm.r2 == Approx(0.8959760903));
    REQUIRE(lm.m  == Approx(0.9451342707));
    REQUIRE(lm.c  == Approx(2.5274566847));

    const auto cmd = "-t FusionExpress -rfus data/fusion/FUS.v1.ref -m data/fusion/FUS.v3.csv -uout tests/data/fusion/fusions.out";
    
    const auto status = parse_options(cmd, _output, _error);
    
    REQUIRE(status == 0);
}