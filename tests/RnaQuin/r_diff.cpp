//#include <catch.hpp>
//#include "unit/test.hpp"
//#include "TransQuin/t_diff.hpp"
//
//using namespace Anaquin;
//
//TEST_CASE("TDiff_DESeq2")
//{
//    Test::transAB();
//
//    TDiff::Options o;
//    
//    o.metrs = TDiff::Metrics::Gene;
//    o.dSoft = TDiff::Software::DESeq2;
//
//    const auto r = TDiff::analyze("tests/data/DESeq2.csv", o);
//
//    REQUIRE(r.size() == 2);
//
//    REQUIRE(r.ids[0]   == "R1_101");
//    REQUIRE(r.ps[0]    == Approx(5.207225e-02));
//    REQUIRE(r.mlfs[0]  == Approx(-1.79218027));
//    REQUIRE(r.elfs[0]  == Approx(-3.0));
//    REQUIRE(r.ses[0]   == Approx(0.92259825));
//    REQUIRE(r.means[0] == Approx(3.924675));
//
//    REQUIRE(r.ids[2]   == "R1_103");
//    REQUIRE(r.ps[2]    == Approx(6.250000e-27));
//    REQUIRE(r.mlfs[2]  == Approx(-0.99634137));
//    REQUIRE(r.elfs[2]  == Approx(-1.0));
//    REQUIRE(r.ses[2]   == Approx(0.09272531));
//    REQUIRE(r.means[2] == Approx(691.727098));
//
//    REQUIRE(r.ids[5]   == "R1_13");
//    REQUIRE(r.ps[5]    == Approx(5.981102e-01));
//    REQUIRE(r.mlfs[5]  == Approx(0.03421518));
//    REQUIRE(r.elfs[5]  == Approx(0.0));
//    REQUIRE(r.ses[5]   == Approx(0.06490966));
//    REQUIRE(r.means[5] == Approx(5658.731489));
//}
//
//TEST_CASE("TDiff_Classify")
//{
//    const auto qvals = std::vector<double> { 0.01, 0.02, 0.98, 0.99 };
//    const auto folds = std::vector<double> { 1.00, 4.00, 1.00, 4.00 };
//
//    const auto r = TDiff::classify(qvals, folds, 0.05, 1.00);
//    
//    REQUIRE(r.size() == 4);
//
//    REQUIRE(r[0] == "FP");
//    REQUIRE(r[1] == "TP");
//    REQUIRE(r[2] == "TN");
//    REQUIRE(r[3] == "FN");
//}
//
//TEST_CASE("TDiff_AllExpressed")
//{
//    const auto &s = Standard::instance().r_trans;
//    
//    Test::transAB();
//    
//    std::set<GeneID> gIDs;
//    
//    for (const auto &i : s.data())
//    {
//        gIDs.insert(i.second.gID);
//    }
//    
//    std::vector<DiffTest> tests;
//
//    for (const auto &gID : gIDs)
//    {
//        DiffTest test;
//        
//        test.id  = gID;
//        test.cID = ChrT;
//        
//        test.p = 0.005;
//        test.q = 0.005;
//
//        // Expected log-fold
//        test.logF = log2(s.findGene(ChrT, gID)->concent(Mix_2) / s.findGene(ChrT, gID)->concent(Mix_1));
//        
//        tests.push_back(test);
//    }
//    
//    TDiff::Options o;
//    
//    o.metrs = TDiff::Metrics::Gene;
//    o.dSoft = TDiff::Software::Cuffdiff;
//    
//    const auto r = TDiff::analyze(tests, o);
//    const auto stats = r.linear();
//    
//    REQUIRE(stats.r  == 1.0);
//    REQUIRE(stats.m  == 1.0);
//    REQUIRE(stats.R2 == 1.0);
//}
//
//TEST_CASE("TDiff_NoneExpressed")
//{
//    const auto &s = Standard::instance().r_trans;
//
//    Test::transAB();
//    
//    std::set<GeneID> gIDs;
//    
//    for (const auto &i : s.data())
//    {
//        gIDs.insert(i.second.gID);
//    }
//    
//    std::vector<DiffTest> tests;
//    
//    for (const auto &gID : gIDs)
//    {
//        DiffTest test;
//        
//        test.id  = gID;
//        test.cID = ChrT;
//        
//        test.p = 0.99;
//        test.q = 0.99;
//        
//        // Expected log-fold
//        test.logF = log2(s.findGene(ChrT, gID)->concent(Mix_2) / s.findGene(ChrT, gID)->concent(Mix_1));
//        
//        tests.push_back(test);
//    }
//    
//    TDiff::Options o;
//    o.metrs = TDiff::Metrics::Gene;
//    
//    const auto r = TDiff::analyze(tests, o);
//    const auto stats = r.linear();
//    
//    REQUIRE(stats.r  == 1.0);
//    REQUIRE(stats.m  == 1.0);
//    REQUIRE(stats.R2 == 1.0);
//}
//
//TEST_CASE("TDiff_NotSynthetic")
//{
//    Test::transAB();
//    
//    std::set<GeneID> gIDs;
//    
//    for (const auto &i : Standard::instance().r_trans.data())
//    {
//        gIDs.insert(i.second.gID);
//    }
//    
//    std::vector<DiffTest> tests;
//    
//    for (const auto &gID : gIDs)
//    {
//        DiffTest test;
//
//        // Randomly generate a unique ID
//        test.id  = gID;
//        
//        // Important, this is telling Anaquin it's not from synthetic chromosome
//        test.cID = "Anaquin";
//        
//        tests.push_back(test);
//    }
//    
//    TDiff::Options o;
//    o.metrs = TDiff::Metrics::Gene;
//    
//    const auto r = TDiff::analyze(tests, o);
//    const auto stats = r.linear();
//    
//    REQUIRE(isnan(stats.r));
//    REQUIRE(isnan(stats.m));
//    REQUIRE(isnan(stats.R2));
//}
