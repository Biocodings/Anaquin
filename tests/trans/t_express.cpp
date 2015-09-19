#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_T_1000_Isoforms")
{
    Test::transA();
    
    TExpress::Options o;
    o.level = Anaquin::TExpress::RNALevel::Isoform;
    
    const auto r  = TExpress::report("tests/data/T_1000/B/G/isoforms.fpkm_tracking", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.r  == Approx(0.232903606));
    REQUIRE(lm.m  == Approx(61.7401661036));
    REQUIRE(lm.r2 == Approx(0.0542440897));
    
    REQUIRE(r.s.id == "R2_38_1");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(0.0000047454));
}

TEST_CASE("TExpress_T_1000_Genes")
{
    Test::transA();
    
    const auto r1 = Test::test("-t TransExpress -m data/trans/MTR002.v013.csv -rgtf data/trans/ATR001.v032.gtf -ugtrack tests/data/T_1000/B/G/genes.fpkm_tracking");
    
    REQUIRE(r1.status == 0);
    
    Test::transA();
    
    TExpress::Options o;
    o.level = Anaquin::TExpress::RNALevel::Gene;
    
    const auto r2 = TExpress::report("tests/data/T_1000/B/G/genes.fpkm_tracking", o);
    const auto lm = r2.linear();
    
    REQUIRE(lm.r  == Approx(0.6079488231));
    REQUIRE(lm.m  == Approx(123.8921186175));
    REQUIRE(lm.r2 == Approx(0.3696017715));
    
    REQUIRE(r2.s.id == "R2_53");
    REQUIRE(r2.s.counts == 1);
    REQUIRE(r2.s.abund == Approx(0.000066038));
}