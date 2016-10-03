#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/r_fold.hpp"

using namespace Anaquin;

TEST_CASE("RnaFoldChange_Sleuth")
{
    Test::RnaQuin_AB();
    
    auto o = RFold::Options();
    
    o.format = RFold::Format::Sleuth;
    o.metrs  = RFold::Metrics::Isoform;
    
    auto r = RFold::analyze("tests/data/sleuth.csv", o);
    
    REQUIRE(r.data.size() == 105);
}

TEST_CASE("RnaFoldChange_Guide")
{
    Test::clear();
    
    const auto r1 = Test::test("RnaFoldChange -o 5.5.3 -m data/RnaQuin/MRN029_v001.csv -method gene -ufiles tests/data/DESeq2.txt");
    const auto r2 = Test::test("RnaFoldChange –o 5.5.3 –m data/RnaQuin/MRN029_v001.csv -method gene -ufiles tests/data/DESeq2.txt");
    
    REQUIRE(r1.status == 0);
    REQUIRE(r2.status == 0);
}
