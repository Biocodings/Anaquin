#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_Guided_Genes")
{
    Test::transA();
 
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Gene;

    auto r = RExpress::analyze("tests/data/guided.gtf", o);

    REQUIRE(r.countSyn == 74);
    REQUIRE(r.countGen == 0);
    REQUIRE(r.dilution() == 1.0);

    REQUIRE(r.gData.size() == 0);
    REQUIRE(r.genes.size() == 74);

    REQUIRE(r.genes["R2_73"].x == Approx(1.8882751465));
    REQUIRE(r.genes["R2_73"].y == Approx(0.3761108));

    const auto ls = r.genes.linear();

    REQUIRE(ls.r  == Approx(0.9520079819));
    REQUIRE(ls.R2 == Approx(0.9063191975));
}

TEST_CASE("TExpress_Guided_Isoforms")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Isoform;
    
    auto r = RExpress::analyze("tests/data/guided.gtf", o);
    
    REQUIRE(r.countSyn == 214);
    REQUIRE(r.countGen == 0);
    REQUIRE(r.dilution() == 1.0);
    
    REQUIRE(r.gData.size() == 0);
    REQUIRE(r.isos.size() == 149);

    REQUIRE(r.isos["R2_73_1"].x == Approx(0.057220459));
    REQUIRE(r.isos["R2_73_1"].y == Approx(0.0000587338));
    REQUIRE(r.isos["R2_73_2"].x == Approx(1.831054688));
    REQUIRE(r.isos["R2_73_2"].y == Approx(0.3760520666));

    const auto ls = r.isos.linear();
    
    REQUIRE(ls.r  == Approx(0.8052543963));
    REQUIRE(ls.R2 == Approx(0.6484346427));
}

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
