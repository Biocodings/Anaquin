#include <catch.hpp>
#include "data/experiment.hpp"

using namespace Anaquin;

TEST_CASE("Experiment_1")
{
    Experiment exp("1,1,1,2,2,2", "A1,A2,A3,B1,B2,B3");

    // The replicates must correspond to the factors in orders
    const auto &samps = exp.samples();

    REQUIRE(samps.size() == 6);
    REQUIRE(samps[0] == 0);
    REQUIRE(samps[1] == 0);
    REQUIRE(samps[2] == 0);
    REQUIRE(samps[3] == 1);
    REQUIRE(samps[4] == 1);
    REQUIRE(samps[5] == 1);
}