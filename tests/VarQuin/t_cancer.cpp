#include <catch.hpp>
#include "test.hpp"

using namespace Anaquin;

TEST_CASE("VarSomatic_1")
{
    REQUIRE_NOTHROW(runTest("VarSomatic -rbed tests/A.V.7.bed -rvcf tests/A.V.8.vcf -useqs tests/sample2.vcf"));
}
