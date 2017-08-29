//#include <catch.hpp>
//#include "test.hpp"
//#include "VarQuin/v_calibrate.hpp"
//
//using namespace Anaquin;
//
//// Defined in resources.cpp
//extern Scripts AVA033Bed();
//
//TEST_CASE("VSubsample_Read_10")
//{
//    clrTest();
//    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
//    
//    VCalibrate::Options o;
//    o.meth  = VCalibrate::Method::Reads;
//    o.reads = 10;
//    
//    auto r = VCalibrate::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
//    
//    REQUIRE(r.count    == 4);
//    
//    REQUIRE(r.beforeEndo == Approx(7.30675));
//    REQUIRE(r.afterEndo  == Approx(7.30675));
//    REQUIRE(r.beforeSeqs == Approx(27.599));
//    REQUIRE(r.afterSeqs  == Approx(7.399));
//    
//    REQUIRE(r.normAver == Approx(0.2576314217));
//    REQUIRE(r.normSD   == Approx(0.2993820454));
//    
//    const auto l1 = Locus(58701157,  58702156);
//    const auto l2 = Locus(60005305,  60006304);
//    const auto l3 = Locus(160204921, 160205920);
//    const auto l4 = Locus(206093670, 206094669);
//    
//    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
//    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.4740740741));
//    REQUIRE(r.c2v["chr1"][l1].endo   == Approx(14.604));
//    REQUIRE(r.c2v["chr1"][l1].before == Approx(29.599));
//    REQUIRE(r.c2v["chr1"][l1].after  == Approx(13.615));
//    
//    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
//    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.5564516129));
//    REQUIRE(r.c2v["chr1"][l2].endo   == Approx(14.623));
//    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.761));
//    REQUIRE(r.c2v["chr1"][l2].after  == Approx(15.981));
//    
//    REQUIRE(r.c2v["chr1"][l3].rID    == "GI_030");
//    REQUIRE(r.c2v["chr1"][l3].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l3].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l3].before == Approx(28.056));
//    REQUIRE(r.c2v["chr1"][l3].after  == 0);
//    
//    REQUIRE(r.c2v["chr1"][l4].rID    == "GI_097");
//    REQUIRE(r.c2v["chr1"][l4].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l4].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l4].before == Approx(25.98));
//    REQUIRE(r.c2v["chr1"][l4].after  == 0);
//}
//
//TEST_CASE("VSubsample_ZeroProp")
//{
//    clrTest();
//    
//    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
//    
//    VCalibrate::Options o;
//    o.meth = VCalibrate::Method::Prop;
//    o.p = 0;
//    
//    auto r = VCalibrate::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
//    
//    REQUIRE(r.count      == 4);
//    REQUIRE(r.afterEndo  == Approx(7.30675));
//    REQUIRE(r.afterSeqs  == 0.0);
//    REQUIRE(r.beforeEndo == Approx(7.30675));
//    REQUIRE(r.beforeSeqs == Approx(27.599));
//
//    REQUIRE(r.normAver == Approx(0.0));
//    REQUIRE(r.normSD   == Approx(0.0));
//    
//    const auto l1 = Locus(58701157,  58702156);
//    const auto l2 = Locus(60005305,  60006304);
//    const auto l3 = Locus(160204921, 160205920);
//    const auto l4 = Locus(206093670, 206094669);
//    
//    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
//    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.0));
//    REQUIRE(r.c2v["chr1"][l1].endo   == Approx(14.604));
//    REQUIRE(r.c2v["chr1"][l1].before == Approx(29.599));
//    REQUIRE(r.c2v["chr1"][l1].after  == Approx(0.0));
//    
//    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
//    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.0));
//    REQUIRE(r.c2v["chr1"][l2].endo   == Approx(14.623));
//    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.761));
//    REQUIRE(r.c2v["chr1"][l2].after  == Approx(0.0));
//    
//    REQUIRE(r.c2v["chr1"][l3].rID    == "GI_030");
//    REQUIRE(r.c2v["chr1"][l3].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l3].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l3].before == Approx(28.056));
//    REQUIRE(r.c2v["chr1"][l3].after  == 0);
//    
//    REQUIRE(r.c2v["chr1"][l4].rID    == "GI_097");
//    REQUIRE(r.c2v["chr1"][l4].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l4].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l4].before == Approx(25.98));
//    REQUIRE(r.c2v["chr1"][l4].after  == 0);
//}
//
//TEST_CASE("VSubsample_Median")
//{
//    clrTest();
//    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
//    
//    VCalibrate::Options o;
//    o.meth = VCalibrate::Method::Median;
//    
//    auto r = VCalibrate::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
//    
//    REQUIRE(r.count      == 4);
//    REQUIRE(r.afterEndo  == Approx(7.0));
//    REQUIRE(r.beforeEndo == Approx(7.0));
//    REQUIRE(r.beforeSeqs == Approx(27.5));
//    
//    REQUIRE(r.afterSeqs >= 7.00);
//    REQUIRE(r.afterSeqs <= 7.70);
//    
//    REQUIRE(r.normAver == Approx(0.2512820513));
//    REQUIRE(r.normSD   == Approx(0.2916321479));
//    
//    const auto l1 = Locus(58701157,  58702156);
//    const auto l2 = Locus(60005305,  60006304);
//    const auto l3 = Locus(160204921, 160205920);
//    const auto l4 = Locus(206093670, 206094669);
//    
//    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
//    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.4666666667));
//    REQUIRE(r.c2v["chr1"][l1].endo   == Approx(14.0));
//    REQUIRE(r.c2v["chr1"][l1].before == Approx(30.0));
//    
//    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
//    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.5384615385));
//    REQUIRE(r.c2v["chr1"][l2].endo   == Approx(14.0));
//    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.0));
//    
//    REQUIRE(r.c2v["chr1"][l3].rID    == "GI_030");
//    REQUIRE(r.c2v["chr1"][l3].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l3].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l3].before == Approx(28.0));
//    REQUIRE(r.c2v["chr1"][l3].after  == 0);
//    
//    REQUIRE(r.c2v["chr1"][l4].rID    == "GI_097");
//    REQUIRE(r.c2v["chr1"][l4].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l4].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l4].before == Approx(26.0));
//    REQUIRE(r.c2v["chr1"][l4].after  == 0);
//}
//
//TEST_CASE("VSubsample_Mean")
//{
//    clrTest();
//    Standard::instance().addVRef(Reader(AVA033Bed(), DataMode::String));
//
//    VCalibrate::Options o;
//    o.meth = VCalibrate::Method::Mean;
//    
//    auto r = VCalibrate::analyze("tests/data/genome.bam", "tests/data/sequins.bam", o);
//    
//    REQUIRE(r.count      == 4);
//    REQUIRE(r.afterEndo  == Approx(7.30675));
//    REQUIRE(r.beforeEndo == Approx(7.30675));
//    REQUIRE(r.beforeSeqs == Approx(27.599));
//
//    REQUIRE(r.afterSeqs >= 7.30);
//    REQUIRE(r.afterSeqs <= 7.70);
//            
//    REQUIRE(r.normAver == Approx(0.2599561382));
//    REQUIRE(r.normSD   == Approx(0.3009513261));
//    
//    const auto l1 = Locus(58701157,  58702156);
//    const auto l2 = Locus(60005305,  60006304);
//    const auto l3 = Locus(160204921, 160205920);
//    const auto l4 = Locus(206093670, 206094669);
//    
//    REQUIRE(r.c2v["chr1"][l1].rID    == "GI_086");
//    REQUIRE(r.c2v["chr1"][l1].norm   == Approx(0.493395));
//    REQUIRE(r.c2v["chr1"][l1].endo   == Approx(14.604));
//    REQUIRE(r.c2v["chr1"][l1].before == Approx(29.599));
//    
//    REQUIRE(r.c2v["chr1"][l2].rID    == "GI_048");
//    REQUIRE(r.c2v["chr1"][l2].norm   == Approx(0.54643));
//    REQUIRE(r.c2v["chr1"][l2].endo   == Approx(14.623));
//    REQUIRE(r.c2v["chr1"][l2].before == Approx(26.761));
//    
//    REQUIRE(r.c2v["chr1"][l3].rID    == "GI_030");
//    REQUIRE(r.c2v["chr1"][l3].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l3].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l3].before == Approx(28.056));
//    REQUIRE(r.c2v["chr1"][l3].after  == 0);
//    
//    REQUIRE(r.c2v["chr1"][l4].rID    == "GI_097");
//    REQUIRE(r.c2v["chr1"][l4].norm   == 0);
//    REQUIRE(r.c2v["chr1"][l4].endo   == 0);
//    REQUIRE(r.c2v["chr1"][l4].before == Approx(25.98));
//    REQUIRE(r.c2v["chr1"][l4].after  == 0);
//}
//
//TEST_CASE("VSubsample_Edge")
//{
//    srand(100);
//    
//    clrTest();
//    Standard::instance().addVRef(Reader("tests/data/test2.bed"));
//    auto r1 = VCalibrate::analyze("tests/data/test3.bam", "tests/data/test2.bam", VCalibrate::Options());
//
//    clrTest();
//    VCalibrate::Options o;
//    o.edge = 300;
//    Standard::instance().addVRef(Reader("tests/data/test2.bed"), o.edge);
//    auto r2 = VCalibrate::analyze("tests/data/test3.bam", "tests/data/test2.bam", o);
//    
//    REQUIRE(r1.beforeEndo == Approx(62.600445186421815));
//    REQUIRE(r1.beforeSeqs == Approx(973.71007234279352));
//    REQUIRE(r1.afterEndo  == Approx(62.600445186421815));
//    REQUIRE(r1.afterSeqs  == Approx(58.5370061213));
//    
//    REQUIRE(r2.beforeEndo == Approx(75.3742690058));
//    REQUIRE(r2.beforeSeqs == Approx(1185.0868838763577));
//    REQUIRE(r2.afterEndo  == Approx(75.3742690058));
//    REQUIRE(r2.afterSeqs  == Approx(60.5175292154));
//}
