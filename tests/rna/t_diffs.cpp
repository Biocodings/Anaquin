#include <catch.hpp>
#include "rna/r_diffs.hpp"

using namespace Spike;

TEST_CASE("Differential_Simulations_Isoforms")
{
    RDiffs::analyze("/Users/tedwong/Sources/QA/scripts/isoform_exp.diff");

    
//    /*
//     * The experiment for mixutre B is identical but it has 100 times coverage.
//     */
//
//    RDiffs::Options o;
//    o.level = RNALevel::Isoform;
//
//    const auto r  = RDiffs::analyze("tests/data/rna/isoform_exp.diff", o);
//    const auto lm = r.linear();
//    
//    REQUIRE(lm.m  == Approx(0.9771063293));
//    REQUIRE(lm.r  == Approx(0.9637073186));
//    REQUIRE(lm.r2 == Approx(0.9268056282));
}

TEST_CASE("Differential_Simulations_Genes")
{
//    /*
//     * The experiment for mixutre B is identical but it has 100 times coverage.
//     */
//
//    RDiffs::Options o;
//    o.level = RNALevel::Gene;
//
//    const auto r  = RDiffs::analyze("tests/data/rna/gene_exp.diff", o);
//    const auto lm = r.linear();
//
//    REQUIRE(lm.m  == Approx(1.0824361534));
//    REQUIRE(lm.r  == Approx(0.8374341202));
//    REQUIRE(lm.r2 == Approx(0.6863607009));
}