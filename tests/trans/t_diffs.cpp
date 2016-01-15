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
        
        test.cID = "chrT";
        test.id  = gID;
        
        test.p = 0.005;
        test.q = 0.005;

        test.logF   = 1.0;
        test.status = DiffTest::Status::NotTested;
        test.fpkm_1 = test.fpkm_2 = 0;
        
        for (const auto &j : Standard::instance().r_trans.data())
        {
            if (gID == j.second.gID)
            {
                test.fpkm_1 += j.second.mixes.at(Mix_1);
                test.fpkm_2 += j.second.mixes.at(Mix_2);
            }
        }
        
        tests.push_back(test);
    }
    
    TDiffs::Options o;
    o.metrs = TDiffs::Metrics::Gene;
    
    const auto r = TDiffs::analyze(tests, o);
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(stats.r  == 1.0);
    REQUIRE(stats.m  == 1.0);
    REQUIRE(stats.r2 == 1.0);
}

TEST_CASE("TDiff_NoneExpressed")
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
        
        test.cID = "chrT";
        test.id  = gID;
        
        test.p = 0.99;
        test.q = 0.99;
        
        test.logF = 1.0;
        test.status = DiffTest::Status::NotTested;
        test.fpkm_1 = test.fpkm_2 = 0;
        
        for (const auto &j : Standard::instance().r_trans.data())
        {
            if (gID == j.second.gID)
            {
                test.fpkm_1 += j.second.mixes.at(Mix_1);
                test.fpkm_2 += j.second.mixes.at(Mix_2);
            }
        }
        
        tests.push_back(test);
    }
    
    TDiffs::Options o;
    o.metrs = TDiffs::Metrics::Gene;
    
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
    o.metrs = TDiffs::Metrics::Gene;
    
    const auto r = TDiffs::analyze(tests, o);
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(isnan(stats.r));
    REQUIRE(isnan(stats.m));
    REQUIRE(isnan(stats.r2));
}