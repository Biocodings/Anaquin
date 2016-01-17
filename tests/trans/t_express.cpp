#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_Perfect")
{
    Test::transA();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        exp.cID  = "chrT";
        exp.id   = i.first;
        exp.fpkm = i.second.mixes.at(Mix_1);
        exps.push_back(exp);
    }

    TExpress::Options o;
    o.lvl = TExpress::Level::Isoform;
    
    const auto r = TExpress::analyze(exps, o);
    const auto stats = r.data.at(ChrT).linear();

    REQUIRE(stats.r  == Approx(1.0));
    REQUIRE(stats.m  == Approx(1.0));
    REQUIRE(stats.r2 == Approx(1.0));
}

TEST_CASE("TExpress_NoSynthetic")
{
    Test::transA();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        // The tool shouldn't process it unless it's chrT
        exp.cID = "Anaquin";
    }
    
    const auto r = TExpress::analyze(exps, TExpress::Options());
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(isnan(stats.r));
    REQUIRE(isnan(stats.m));
    REQUIRE(isnan(stats.r2));
}

TEST_CASE("TExpress_AllZeros")
{
    Test::transA();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        exp.cID  = "chrT";
        exp.id   = i.first;
        exp.fpkm = 0.0;
        exps.push_back(exp);
    }
    
    const auto r = TExpress::analyze(exps, TExpress::Options());
    const auto stats = r.data.at(ChrT).linear();
    
    REQUIRE(isnan(stats.r));
    REQUIRE(isnan(stats.m));
    REQUIRE(isnan(stats.r2));
}

TEST_CASE("TExpress_FlatMix")
{
    Test::transF();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        exp.cID  = "chrT";
        exp.id   = i.first;
        exp.fpkm = 100 * rand();
        exps.push_back(exp);
    }
    
    const auto r = TExpress::analyze(exps, TExpress::Options());
    const auto stats = r.data.at(ChrT).linear();

    REQUIRE(isnan(stats.r));
    REQUIRE(isnan(stats.m));
    REQUIRE(isnan(stats.r2));
}
