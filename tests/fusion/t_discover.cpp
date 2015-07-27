#include <catch.hpp>
#include "fusion/f_discover.hpp"

static std::string _output, _error;

extern int parse_options(const std::string &, std::string &, std::string &);

using namespace Anaquin;

TEST_CASE("FDiscover_Command")
{
    const auto cmd = "-t FusionDiscover -rfus data/fusion/FUS.v1.ref -m data/fusion/FUS.v3.csv -uout tests/data/fusion/fusions.out";
    const auto status = parse_options(cmd, _output, _error);

    REQUIRE(status == 0);
}

TEST_CASE("FDiscover_Test")
{
    const auto stats = FDiscover::analyze("tests/data/fusion/fusions.out");

    REQUIRE(stats.m.sn() == Approx(0.7916666667));
    REQUIRE(stats.m.sp() == Approx(1.0));
    REQUIRE(stats.m.nr == 24);
}