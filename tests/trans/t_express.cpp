#include <catch.hpp>
#include "trans/t_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_Simulations_Genes_A1")
{
    TExpress::Options opt;
    opt.level = Anaquin::TExpress::RNALevel::Gene;

    const auto r  = TExpress::analyze("tests/data/rna/A1/genes.fpkm_tracking", opt);
    const auto lm = r.linear();

    REQUIRE(lm.r  == Approx(0.8929552651));
    REQUIRE(lm.m  == Approx(0.9065369839));
    REQUIRE(lm.r2 == Approx(0.7920367134));
    
    REQUIRE(r.s.id == "R2_33");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(7.0));
}

TEST_CASE("TExpress_Simulations_Isoforms_A1")
{
    const auto r  = TExpress::analyze("tests/data/rna/A1/isoforms.fpkm_tracking");
    const auto lm = r.linear();

    REQUIRE(lm.r  == Approx(0.8143944753));
    REQUIRE(lm.m  == Approx(0.9300431127));
    REQUIRE(lm.r2 == Approx(0.6591315122));

    REQUIRE(r.s.id == "R2_38_1");
    REQUIRE(r.s.counts == 1);
    REQUIRE(r.s.abund == Approx(1.0));
}