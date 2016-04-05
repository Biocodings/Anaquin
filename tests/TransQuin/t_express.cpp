#include <catch.hpp>
#include "unit/test.hpp"
#include "TransQuin/t_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_AllZeros")
{
    Test::transA();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        exp.cID   = "chrT";
        exp.id    = i.first;
        exp.abund = 0.0;
        exps.push_back(exp);
    }
    
    TExpress::Options o;
    o.metrs = TExpress::Metrics::Isoform;
    
    REQUIRE_THROWS(TExpress::analyze(exps, o));
}

TEST_CASE("TExpress_Perfect")
{
    Test::transA();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        exp.cID   = "chrT";
        exp.id    = i.first;
        exp.abund = i.second.mixes.at(Mix_1);
        exps.push_back(exp);
    }

    TExpress::Options o;
    o.metrs = TExpress::Metrics::Isoform;
    
    const auto r = TExpress::analyze(exps, o);
    const auto stats = r.data.linear();

    REQUIRE(stats.r  == Approx(1.0));
    REQUIRE(stats.m  == Approx(1.0));
    REQUIRE(stats.R2 == Approx(1.0));
}

TEST_CASE("TExpress_NoSynthetic")
{
    Test::transA();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;        
        exp.cID = "Anaquin";
    }
    
    REQUIRE_THROWS(TExpress::analyze(exps, TExpress::Options()));
}

TEST_CASE("TExpress_FlatMix")
{
    Test::transF();
    
    std::vector<Expression> exps;
    
    for (const auto &i : Standard::instance().r_trans.data())
    {
        Expression exp;
        
        exp.cID   = "chrT";
        exp.id    = i.first;
        exp.abund = 100 * rand();
        exps.push_back(exp);
    }
    
    TExpress::Options o;
    o.metrs = TExpress::Metrics::Isoform;
    
    const auto r = TExpress::analyze(exps, o);
    const auto stats = r.data.linear();

    REQUIRE(isnan(stats.r));
    REQUIRE(isnan(stats.m));
    REQUIRE(isnan(stats.R2));
}
