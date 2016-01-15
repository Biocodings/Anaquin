#include <catch.hpp>
#include "data/experiment.hpp"

using namespace Anaquin;

TEST_CASE("Experiment_1")
{
    Experiment exp("1,1,1,2,2,2");

    // The replicates must correspond to the factors in orders
    const auto &reps = exp.reps();

    REQUIRE(reps.size() == 6);
    REQUIRE(reps[0] == 0);
    REQUIRE(reps[1] == 0);
    REQUIRE(reps[2] == 0);
    REQUIRE(reps[3] == 1);
    REQUIRE(reps[4] == 1);
    REQUIRE(reps[5] == 1);
}