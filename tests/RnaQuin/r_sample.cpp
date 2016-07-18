#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_sample.hpp"

using namespace Anaquin;

TEST_CASE("RSample_Negative")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method -0.5 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. -0.5 not expected for option -method. Sampling fraction must be greater than zero.\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Zero")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 0.00 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. method. Sampling fraction must be greater than zero not expected for -0.00.\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_One")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 1.00 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. method. Sampling fraction must be less than one not expected for -1.00.\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Ten")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 10 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "Invalid command. 10 not expected for option -method. Sampling fraction must be less than one.\n");
    REQUIRE(r.status == 1);
}