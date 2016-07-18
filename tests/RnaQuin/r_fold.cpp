#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_fold.hpp"

using namespace Anaquin;

TEST_CASE("RnaFoldChange_Manuscript")
{
    Test::clear();
    
    const auto r = Test::test("RnaFoldChange -o 5.5.3 -m data/RnaQuin/MRN029_v001.csv -method gene -ufiles tests/data/DESeq2.txt");

    REQUIRE(r.status == 0);
}