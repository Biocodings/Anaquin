#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_diffs.hpp"

using namespace Anaquin;

TEST_CASE("TDiff_Classify")
{
    const auto qvals = std::vector<double> { 0.01, 0.02, 0.98, 0.99 };
    const auto folds = std::vector<double> { 1.00, 4.00, 1.00, 4.00 };

    const auto r = TDiffs::classify(qvals, folds, 0.05, 1.00);
    
    REQUIRE(r.size() == 4);

    REQUIRE(r[0] == "FP");
    REQUIRE(r[1] == "TP");
    REQUIRE(r[2] == "TN");
    REQUIRE(r[3] == "FN");
}

TEST_CASE("TDiff_AllExpressed")
{
    const auto &s = Standard::instance().r_trans;
    
    Test::transAB();
    
    std::set<GeneID> gIDs;
    
    for (const auto &i : s.data())
    {
        gIDs.insert(i.second.gID);
    }
    
    std::vector<DiffTest> tests;

    for (const auto &gID : gIDs)
    {
        DiffTest test;
        
        test.id  = gID;
        test.cID = ChrT;
        
        test.p = 0.005;
        test.q = 0.005;

        // Expected log-fold
        test.logF = log2(s.findGene(ChrT, gID)->abund(Mix_2) / s.findGene(ChrT, gID)->abund(Mix_1));
        
        tests.push_back(test);
    }
    
    TDiffs::Options o;
    
    o.lvl   = TDiffs::Level::Gene;
    o.dSoft = TDiffs::DiffSoft::Cuffdiff;
    
    const auto r = TDiffs::analyze(tests, o);
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(stats.r  == 1.0);
    REQUIRE(stats.m  == 1.0);
    REQUIRE(stats.r2 == 1.0);
}

TEST_CASE("TDiff_NoneExpressed")
{
    const auto &s = Standard::instance().r_trans;

    Test::transAB();
    
    std::set<GeneID> gIDs;
    
    for (const auto &i : s.data())
    {
        gIDs.insert(i.second.gID);
    }
    
    std::vector<DiffTest> tests;
    
    for (const auto &gID : gIDs)
    {
        DiffTest test;
        
        test.id  = gID;
        test.cID = ChrT;
        
        test.p = 0.99;
        test.q = 0.99;
        
        // Expected log-fold
        test.logF = log2(s.findGene(ChrT, gID)->abund(Mix_2) / s.findGene(ChrT, gID)->abund(Mix_1));
        
        tests.push_back(test);
    }
    
    TDiffs::Options o;
    o.lvl = TDiffs::Level::Gene;
    
    const auto r = TDiffs::analyze(tests, o);
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(stats.r  == 1.0);
    REQUIRE(stats.m  == 1.0);
    REQUIRE(stats.r2 == 1.0);
}

TEST_CASE("TDiff_NoSynthetic")
{
    Test::transAB();
    
    std::set<GeneID> gIDs;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        gIDs.insert(i.second.gID);
    }
    
    std::vector<DiffTest> tests;
    
    for (const auto &gID : gIDs)
    {
        DiffTest test;
        test.cID = "Anaquin";
        tests.push_back(test);
    }
    
    TDiffs::Options o;
    o.lvl = TDiffs::Level::Gene;
    
    const auto r = TDiffs::analyze(tests, o);
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(isnan(stats.r));
    REQUIRE(isnan(stats.m));
    REQUIRE(isnan(stats.r2));
}