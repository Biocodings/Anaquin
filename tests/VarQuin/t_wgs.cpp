#include <catch.hpp>
#include "test.hpp"

using namespace Anaquin;

TEST_CASE("VarWGS_1")
{
    REQUIRE_NOTHROW(runTest("VarWGS -rvcf tests/A.V.6.vcf -rbed tests/A.V.5.bed -useqs tests/sample1.vcf"));
}
