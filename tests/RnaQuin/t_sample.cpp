#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/r_sample.hpp"

using namespace Anaquin;

TEST_CASE("RSample_OverSubsampled")
{
    clrTest();    
    RSample::Options o;

    // Over subsampling
    o.p = 0.99;
    
    const auto r = RSample::stats("tests/data/sampled.bam", o);
    
    REQUIRE(r.before.syn == 32242);
    REQUIRE(r.before.gen == 420874);
    REQUIRE(r.before.dilut() == Approx(0.0711561719));
    
    REQUIRE(r.after.syn == 32242);
    REQUIRE(r.after.gen == 420874);
    REQUIRE(r.after.dilut() == Approx(0.0711561719));
}

TEST_CASE("RSample_Negative")
{
    clrTest();
    const auto r = runTest("RnaSubsample -method -0.5 -usequin tests/data/sampled.bam");
    
    REQUIRE(r.error == "***********************\n[ERRO]: Invalid value for -method. Sampling fraction must be greater than zero.\n***********************\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Zero")
{
    clrTest();
    const auto r = runTest("RnaSubsample -method 0.00 -usequin tests/data/sampled.bam");
    
    REQUIRE(r.error == "***********************\n[ERRO]: Invalid value for -method. Sampling fraction must be greater than zero.\n***********************\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_One")
{
    clrTest();
    const auto r = runTest("RnaSubsample -method 1.00 -usequin tests/data/sampled.bam");
    
    REQUIRE(r.error == "***********************\n[ERRO]: Invalid value for -method. Sampling fraction must be less than one.\n***********************\n");
    REQUIRE(r.status == 1);
}

TEST_CASE("RSample_Ten")
{
    clrTest();
    const auto r = runTest("RnaSubsample -method 10 -usequin tests/data/sampled.bam");
    
    REQUIRE(r.error == "***********************\n[ERRO]: Invalid value for -method. Sampling fraction must be less than one.\n***********************\n");
    REQUIRE(r.status == 1);
}

