#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/r_fold.hpp"

using namespace Anaquin;

TEST_CASE("RnaFoldChange_Sleuth_1")
{
    RnaQuin_AB();
    
    auto o = RFold::Options();
    
    o.format = RFold::Format::Sleuth;
    o.metrs  = RFold::Metrics::Gene;
    
    auto r = RFold::analyze("tests/data/sleuth.csv", o);
    
    REQUIRE(r.data.size()  == 60);
    REQUIRE(r.means.size() == 0);
    REQUIRE(r.ses.size()   == 0);
    
    REQUIRE(r.data.count("R2_60"));
    REQUIRE(r.data["R2_60"].exp == 4.0);
    REQUIRE(r.data["R2_60"].obs == Approx(3.721217));
    REQUIRE(isnan(r.data["R2_60"].samp1));
    REQUIRE(isnan(r.data["R2_60"].samp2));
    REQUIRE(isnan(r.data["R2_60"].se));
    REQUIRE(isnan(r.data["R2_60"].mean));
    REQUIRE(isnan(r.data["R2_60"].p));
    REQUIRE(isnan(r.data["R2_60"].q));
}

TEST_CASE("RnaFoldChange_Sleuth_2")
{
    RnaQuin_AB();
    
    auto o = RFold::Options();
    
    o.format = RFold::Format::Sleuth;
    o.metrs  = RFold::Metrics::Isoform;
    
    auto r = RFold::analyze("tests/data/sleuth.csv", o);
    
    REQUIRE(r.data.size()  == 105);
    REQUIRE(r.means.size() == 0);
    REQUIRE(r.ses.size()   == 0);
    
    REQUIRE(r.data["R2_60_1"].exp == 8.0);
    REQUIRE(r.data["R2_60_1"].obs == Approx(3.85074248350044));
    REQUIRE(isnan(r.data["R2_60_1"].samp1));
    REQUIRE(isnan(r.data["R2_60_1"].samp2));
    REQUIRE(r.data["R2_60_1"].se == Approx(0.522653700140967));
    REQUIRE(r.data["R2_60_1"].mean == Approx(4.85372795921658));
    REQUIRE(r.data["R2_60_1"].p == Approx(1.73629727192869e-13));
    REQUIRE(r.data["R2_60_1"].q == Approx(4.34074317982173e-13));
}

TEST_CASE("RnaFoldChange_Guide")
{
    clrTest();
    
    const auto r1 = runTest("RnaFoldChange -o 5.5.3 -m data/RnaQuin/MRN029_v001.csv -method gene -ufiles tests/data/DESeq2.txt");
    const auto r2 = runTest("RnaFoldChange –o 5.5.3 –m data/RnaQuin/MRN029_v001.csv -method gene -ufiles tests/data/DESeq2.txt");
    
    REQUIRE(r1.status == 0);
    REQUIRE(r2.status == 0);
}
