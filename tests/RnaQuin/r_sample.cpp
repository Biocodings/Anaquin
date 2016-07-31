#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_sample.hpp"

using namespace Anaquin;

TEST_CASE("RSample_OverSubsampled")
{
    Test::clear();
    
    RSample::Options o;

    // Over subsampling
    o.p = 0.99;
    
    // It's a unit test...
    o.toConsole = false;
    
    const auto r = RSample::stats("tests/data/sampled.bam", o);
    
    REQUIRE(r.before.syn == 32323);
    REQUIRE(r.before.gen == 424185);
    REQUIRE(r.before.dilut() == Approx(0.0708048928));
    
    REQUIRE(r.after.syn == 32323);
    REQUIRE(r.after.gen == 424185);
    REQUIRE(r.after.dilut() == Approx(0.0708048928));
}

TEST_CASE("RSample_Negative")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method -0.5 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "[ERRO]: Invalid value for -method. Sampling fraction must be greater than zero.\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Zero")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 0.00 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "[ERRO]: Invalid value for -method. Sampling fraction must be greater than zero.\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_One")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 1.00 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "[ERRO]: Invalid value for -method. Sampling fraction must be less than one.\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Ten")
{
    Test::clear();
    
    const auto r = Test::test("RnaSubsample -method 10 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.error == "[ERRO]: Invalid value for -method. Sampling fraction must be less than one.\n");
    REQUIRE(r.status == 1);
}