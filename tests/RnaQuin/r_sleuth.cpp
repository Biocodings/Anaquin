#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_fold.hpp"

using namespace Anaquin;

TEST_CASE("RnaSleuth_Workflow")
{
    Test::clear();
    
    const auto r = Test::test("RnaFoldChange -o 5.5.3 -m data/RnaQuin/MTR004.v013.csv -ufiles tests/data/sleuth.txt");
    
    REQUIRE(r.status == 0);
}
