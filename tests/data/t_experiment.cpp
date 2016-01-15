#include <catch.hpp>
#include "data/experiment.hpp"

using namespace Anaquin;

TEST_CASE("Experiment_1")
{
    Experiment exp("1,1,1,2,2,2");

    const auto &facts = exp.factors();
    
    REQUIRE(facts.size() == 6);
    REQUIRE(facts[0] == 0);
    REQUIRE(facts[1] == 0);
    REQUIRE(facts[2] == 0);
    REQUIRE(facts[3] == 1);
    REQUIRE(facts[4] == 1);
    REQUIRE(facts[5] == 1);
}