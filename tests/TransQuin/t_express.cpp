//#include <catch.hpp>
//#include "unit/test.hpp"
//#include "TransQuin/t_express.hpp"
//
//using namespace Anaquin;
//
//TEST_CASE("TExpress_AllZeros")
//{
//    Test::transA();
//    
//    std::vector<ParserCufflink::Data> x;
//    
//    for (const auto &i : Standard::instance().r_trans.data())
//    {
//        TExpress::TestData d;
//        
//        d.cID   = "chrT";
//        d.id    = i.first;
//        d.abund = 0.0;
//        x.push_back(d);
//    }
//    
//    TExpress::Options o;
//    o.metrs = TExpress::Metrics::Isoform;
//    
//    REQUIRE_THROWS(TExpress::analyze(x, o));
//}
//
//TEST_CASE("TExpress_Perfect")
//{
//    Test::transA();
//    
//    std::vector<ParserCufflink::Data> x;
//    
//    for (const auto &i : Standard::instance().r_trans.data())
//    {
//        TExpress::TestData d;
//        
//        d.cID   = "chrT";
//        d.id    = i.first;
//        d.abund = i.second.mixes.at(Mix_1);
//        x.push_back(d);
//    }
//
//    TExpress::Options o;
//    o.metrs = TExpress::Metrics::Isoform;
//    
//    const auto r = TExpress::analyze(x, o);
//    const auto stats = r.data.linear();
//
//    REQUIRE(stats.r  == Approx(1.0));
//    REQUIRE(stats.m  == Approx(1.0));
//    REQUIRE(stats.R2 == Approx(1.0));
//}
//
//TEST_CASE("TExpress_NoSynthetic")
//{
//    Test::transA();
//    
//    std::vector<TExpress::TestData> exps;
//    
//    for (const auto &i : Standard::instance().r_trans.data())
//    {
//        TExpress::TestData d;
//        d.cID = "Anaquin";
//    }
//    
//    REQUIRE_THROWS(TExpress::analyze(exps, TExpress::Options()));
//}
//
//TEST_CASE("TExpress_FlatMix")
//{
//    Test::transF();
//    
//    std::vector<TExpress::TestData> x;
//    
//    for (const auto &i : Standard::instance().r_trans.data())
//    {
//        TExpress::TestData d;
//        
//        d.cID   = "chrT";
//        d.id    = i.first;
//        d.abund = 100 * rand();
//        x.push_back(d);
//    }
//    
//    TExpress::Options o;
//    o.metrs = TExpress::Metrics::Isoform;
//    
//    const auto r = TExpress::analyze(x, o);
//    const auto stats = r.data.linear();
//
//    REQUIRE(isnan(stats.r));
//    REQUIRE(isnan(stats.m));
//    REQUIRE(isnan(stats.R2));
//}
