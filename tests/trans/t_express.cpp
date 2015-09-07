#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_T_1000_Genes")
{
    Test::trans();

    const auto r1 = Test::test("-t TransExpress -m data/trans/TransMixture_4.1.csv -rgtf data/trans/TransStandard_1.0.gtf -ugtrack tests/data/T_1000/B/G/genes.fpkm_tracking");

    REQUIRE(r1.status == 0);

    Test::trans();

    TExpress::Options o;
    o.level = Anaquin::TExpress::RNALevel::Gene;

    const auto r2 = TExpress::analyze("tests/data/T_1000/B/G/genes.fpkm_tracking", o);
    const auto lm = r2.linear();

    REQUIRE(lm.r  == Approx(0.6294323776));
    REQUIRE(lm.m  == Approx(0.592652157));
    REQUIRE(lm.r2 == Approx(0.3836056413));
    
    REQUIRE(r2.s.id == "R2_33");
    REQUIRE(r2.s.counts == 1);
    REQUIRE(r2.s.abund == Approx(7.0));
}

TEST_CASE("TExpress_T_1000_Isoforms")
{
    TExpress::Options o;
    o.level = Anaquin::TExpress::RNALevel::Isoform;

    const auto r  = TExpress::analyze("tests/data/T_1000/B/G/genes.fpkm_tracking", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.r  == Approx(0.8929552651));
    REQUIRE(lm.m  == Approx(0.9065369839));
    REQUIRE(lm.r2 == Approx(0.7920367134));

    REQUIRE(r.s.id == "R2_33");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(7.0));
}