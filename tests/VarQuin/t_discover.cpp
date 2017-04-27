#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

TEST_CASE("VDiscover_1")
{

}

TEST_CASE("VDiscover_2")
{
    clrTest();
    
    REQUIRE_NOTHROW(runTest("VarDiscover -rbed tests/A.V.6.bed -rvcf tests/A.V.8.vcf -ufiles tests/sample.vcf -ufiles tests/sample.vcf"));
}
