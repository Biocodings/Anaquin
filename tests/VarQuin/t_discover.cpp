#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_align.hpp"

using namespace Anaquin;

TEST_CASE("VDiscover_ToyExample")
{
    Test::clear();
    
    REQUIRE_NOTHROW(Test::test("VarDiscover -m tests/M.R.8.csv -rbed tests/A.V.5.bed -rvcf tests/A.V.8.vcf -ufiles tests/sample.vcf"));
}
