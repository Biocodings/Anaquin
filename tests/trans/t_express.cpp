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
    
    REQUIRE(lm.r  == Approx(0.5303389953));
    REQUIRE(lm.m  == Approx(0.4560194933));
    REQUIRE(lm.r2 == Approx(0.2812594499));
    
    REQUIRE(r.ss.id == "R2_38_1");
    REQUIRE(r.ss.counts == 1);
    REQUIRE(r.ss.abund == Approx(0.003933907));
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
    
    REQUIRE(lm.r  == Approx(0.6447104779));
    REQUIRE(lm.m  == Approx(0.6184389099));
    REQUIRE(lm.r2 == Approx(0.4156516003));
    
    REQUIRE(r2.ss.id == "R2_33");
    REQUIRE(r2.ss.counts == 1);
    REQUIRE(r2.ss.abund == Approx(0.0590085983));
}