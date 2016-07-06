#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_sample.hpp"

using namespace Anaquin;

TEST_CASE("RSample_Negative")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method -0.5 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. -0.5 not expected for option -method. Sampling fraction must be greater than zero.");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Zero")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 0.00 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. Unknown tool: RnaSample. Please check the user manual and try again.");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_One")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 1.00 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. Unknown tool: RnaSample. Please check the user manual and try again.");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Ten")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 10 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. Unknown tool: RnaSample. Please check the user manual and try again.");
    REQUIRE(r.status == 1);
}