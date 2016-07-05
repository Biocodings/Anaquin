#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_sample.hpp"

using namespace Anaquin;

TEST_CASE("RSample_QuickCheck")
{
    Test::clear();
    
    const auto r = Test::test("RnaSample -method -25 -ufiles tests/data/sampled.bam");
    
    REQUIRE(r.status == 0);
}