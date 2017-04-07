#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_sample.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts AVA033Bed();

TEST_CASE("VSubsample_Read_10")
{
    Test::clear();
    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
    
    VSample::Options o;
    o.meth  = VSample::Method::Reads;
    o.reads = 10;
    
    auto r = VSample::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
    
    REQUIRE(r.count    == 4);
    REQUIRE(r.noGAlign == 2);
    REQUIRE(r.noSAlign == 0);
    
    REQUIRE(r.beforeGen == Approx(7.30675));
    REQUIRE(r.afterGen  == Approx(7.30675));
    REQUIRE(r.beforeSyn == Approx(27.599));
    REQUIRE(r.afterSyn  == Approx(7.399));
    
    REQUIRE(r.normAver == Approx(0.2576314217));
    REQUIRE(r.normSD   == Approx(0.2993820454));
    
    const auto l1 = Locus(58701157,  58702156);
    const auto l2 = Locus(60005305,  60006304);
    const auto l3 = Locus(160204921, 160205920);
    const auto l4 = Locus(206093670, 206094669);
    
    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.4740740741));
    REQUIRE(r.c2v["chr1"][l1].gen    == Approx(14.604));
    REQUIRE(r.c2v["chr1"][l1].before == Approx(29.599));
    REQUIRE(r.c2v["chr1"][l1].after  == Approx(13.615));
    
    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.5564516129));
    REQUIRE(r.c2v["chr1"][l2].gen    == Approx(14.623));
    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.761));
    REQUIRE(r.c2v["chr1"][l2].after  == Approx(15.981));
    
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

TEST_CASE("VSubsample_ZeroProp")
{
    Test::clear();
    
    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
    
    VSample::Options o;
    o.meth = VSample::Method::Prop;
    o.p = 0;
    
    auto r = VSample::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
    
    REQUIRE(r.count    == 4);
    REQUIRE(r.noGAlign == 2);
    REQUIRE(r.noSAlign == 0);
    
    REQUIRE(r.afterGen  == Approx(7.30675));
    REQUIRE(r.afterSyn  == 0.0);
    REQUIRE(r.beforeGen == Approx(7.30675));
    REQUIRE(r.beforeSyn == Approx(27.599));

    REQUIRE(r.normAver == Approx(0.0));
    REQUIRE(r.normSD   == Approx(0.0));
    
    const auto l1 = Locus(58701157,  58702156);
    const auto l2 = Locus(60005305,  60006304);
    const auto l3 = Locus(160204921, 160205920);
    const auto l4 = Locus(206093670, 206094669);
    
    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.0));
    REQUIRE(r.c2v["chr1"][l1].gen    == Approx(14.604));
    REQUIRE(r.c2v["chr1"][l1].before == Approx(29.599));
    REQUIRE(r.c2v["chr1"][l1].after  == Approx(0.0));
    
    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.0));
    REQUIRE(r.c2v["chr1"][l2].gen    == Approx(14.623));
    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.761));
    REQUIRE(r.c2v["chr1"][l2].after  == Approx(0.0));
    
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

TEST_CASE("VSubsample_Median")
{
    Test::clear();
    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
    
    VSample::Options o;
    o.meth = VSample::Method::Median;
    
    auto r = VSample::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
    
    REQUIRE(r.count    == 4);
    REQUIRE(r.noGAlign == 2);
    REQUIRE(r.noSAlign == 0);
    
    REQUIRE(r.afterGen  == Approx(7.0));
    REQUIRE(r.beforeGen == Approx(7.0));
    REQUIRE(r.beforeSyn == Approx(27.5));
    
    REQUIRE(r.afterSyn >= 7.00);
    REQUIRE(r.afterSyn <= 7.70);
    
    REQUIRE(r.normAver == Approx(0.2512820513));
    REQUIRE(r.normSD   == Approx(0.2916321479));
    
    const auto l1 = Locus(58701157,  58702156);
    const auto l2 = Locus(60005305,  60006304);
    const auto l3 = Locus(160204921, 160205920);
    const auto l4 = Locus(206093670, 206094669);
    
    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.4666666667));
    REQUIRE(r.c2v["chr1"][l1].gen    == Approx(14.0));
    REQUIRE(r.c2v["chr1"][l1].before == Approx(30.0));
    
    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.5384615385));
    REQUIRE(r.c2v["chr1"][l2].gen    == Approx(14.0));
    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.0));
    
    REQUIRE(r.c2v["chr1"][l3].rID    == "GI_030");
    REQUIRE(r.c2v["chr1"][l3].norm   == 0);
    REQUIRE(r.c2v["chr1"][l3].gen    == 0);
    REQUIRE(r.c2v["chr1"][l3].before == Approx(28.0));
    REQUIRE(r.c2v["chr1"][l3].after  == 0);
    
    REQUIRE(r.c2v["chr1"][l4].rID    == "GI_097");
    REQUIRE(r.c2v["chr1"][l4].norm   == 0);
    REQUIRE(r.c2v["chr1"][l4].gen    == 0);
    REQUIRE(r.c2v["chr1"][l4].before == Approx(26.0));
    REQUIRE(r.c2v["chr1"][l4].after  == 0);
}

TEST_CASE("VSubsample_Mean")
{
    Test::clear();
    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));

    VSample::Options o;
    o.meth = VSample::Method::Mean;
    
    auto r = VSample::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
    
    REQUIRE(r.count    == 4);
    REQUIRE(r.noGAlign == 2);
    REQUIRE(r.noSAlign == 0);
    
    REQUIRE(r.afterGen  == Approx(7.30675));
    REQUIRE(r.beforeGen == Approx(7.30675));
    REQUIRE(r.beforeSyn == Approx(27.599));

    REQUIRE(r.afterSyn >= 7.30);
    REQUIRE(r.afterSyn <= 7.70);
            
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

TEST_CASE("VSubsample_Edge")
{
    srand(100);
    
    Test::clear();
    Standard::instance().addVRef(Reader("tests/data/test2.bed"));
    auto r1 = VSample::analyze("tests/data/test3.bam", "tests/data/test2.bam", VSample::Options());

    Test::clear();
    VSample::Options o;
    o.edge = 300;
    Standard::instance().addVRef(Reader("tests/data/test2.bed"), o.edge);
    auto r2 = VSample::analyze("tests/data/test3.bam", "tests/data/test2.bam", o);
    
    REQUIRE(r1.beforeGen == Approx(62.600445186421815));
    REQUIRE(r1.beforeSyn == Approx(973.71007234279352));
    REQUIRE(r1.afterGen  == Approx(62.600445186421815));
    REQUIRE(r1.afterSyn  == Approx(58.5370061213));
    
    REQUIRE(r2.beforeGen == Approx(75.3742690058));
    REQUIRE(r2.beforeSyn == Approx(1185.0868838763577));
    REQUIRE(r2.afterGen  == Approx(75.3742690058));
    REQUIRE(r2.afterSyn  == Approx(60.5175292154));
}
