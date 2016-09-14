#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_sample2.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts AVA033Bed();

// Defined in main.cpp
extern bool __hackBedFile__;

TEST_CASE("VSubsample_Mean")
{
    __hackBedFile__ = true;
    
    Test::clear();

    Standard::instance().addVStd(Reader(AVA033Bed(), DataMode::String));

    VSample2::Options o;
    o.meth = VSample2::Method::Mean;
    
    auto r = VSample2::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
    
    REQUIRE(r.count == 4);
    REQUIRE(r.noGAlign == 2);
    REQUIRE(r.noSAlign == 0);
    
    REQUIRE(r.afterGen  == Approx(7.30675));
    REQUIRE(r.beforeGen == Approx(7.30675));
    REQUIRE(r.beforeSyn == Approx(27.599));

    REQUIRE(r.afterSyn >= 7.36);
    REQUIRE(r.afterSyn <= 7.37);
            
    REQUIRE(r.normAver == Approx(0.2599561382));
    REQUIRE(r.normSD   == Approx(0.3009513261));
    
    const auto l1 = Locus(58701157,  58702156);
    const auto l2 = Locus(60005305,  60006304);
    const auto l3 = Locus(160204921, 160205920);
    const auto l4 = Locus(206093670, 206094669);
    
    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.493395));
    REQUIRE(r.c2v["chr1"][l1].gen    == Approx(14.604));
    REQUIRE(r.c2v["chr1"][l1].before == Approx(29.599));
    
    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.54643));
    REQUIRE(r.c2v["chr1"][l2].gen    == Approx(14.623));
    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.761));
    
    REQUIRE(r.c2v["chr1"][l3].rID    == "GI_030");
    REQUIRE(r.c2v["chr1"][l3].norm   == 0);
    REQUIRE(r.c2v["chr1"][l3].gen    == 0);
    REQUIRE(r.c2v["chr1"][l3].before == Approx(28.056));
    REQUIRE(r.c2v["chr1"][l3].after  == 0);
    
    REQUIRE(r.c2v["chr1"][l4].rID    == "GI_097");
    REQUIRE(r.c2v["chr1"][l4].norm   == 0);
    REQUIRE(r.c2v["chr1"][l4].gen    == 0);
    REQUIRE(r.c2v["chr1"][l4].before == Approx(25.98));
    REQUIRE(r.c2v["chr1"][l4].after  == 0);
}