#include <catch.hpp>
#include "trans/t_diffs.hpp"

using namespace Anaquin;

TEST_CASE("TDiffs_T_1000_Genes")
{
    TDiffs::Options o;
    o.level = TDiffs::Gene;

    const auto r  = TDiffs::analyze("tests/data/T_1000/gene_exp.diff", o);
//    const auto lm = r.linear();
//    
//    REQUIRE(lm.m  == Approx(0.9771063293));
//    REQUIRE(lm.r  == Approx(0.9637073186));
//    REQUIRE(lm.r2 == Approx(0.9268056282));
}

TEST_CASE("TDiffs_T_1000_Isoforms")
{
    TDiffs::Options o;
    o.level = TDiffs::Isoform;
    
    const auto r  = TDiffs::analyze("tests/data/T_1000/isoform_exp.diff", o);
    const auto lm = r.linear();

    REQUIRE(lm.m  == Approx(1.0824361534));
    REQUIRE(lm.r  == Approx(0.8374341202));
    REQUIRE(lm.r2 == Approx(0.6863607009));
}