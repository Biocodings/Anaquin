#include <catch.hpp>
#include "rna/r_differential.hpp"

using namespace Spike;

TEST_CASE("Differential_Simulations_Isoforms")
{
    /*
     * The experiment for mixutre B is identical but it has 20 times coverage.
     */
    
    RDifferential::Options o;
    o.level = RNALevel::Isoform;
    auto r = RDifferential::analyze("tests/data/rna_sims/isoform_exp.diff", o);

    REQUIRE(r.lm.m  == Approx(1.24));
    REQUIRE(r.lm.r  == Approx(1.24));
    REQUIRE(r.lm.r2 == Approx(1.24));
    REQUIRE(r.lm.r2 == Approx(1.24));
}

TEST_CASE("Differential_Simulations_Genes")
{
    /*
     * The experiment for mixutre B is identical but it has 20 times coverage.
     */
    
    RDifferential::Options o;
    o.level = RNALevel::Gene;
    const auto r = RDifferential::analyze("tests/data/rna_sims/gene_exp.diff", o);

    REQUIRE(r.lm.m  == Approx(1.24));
    REQUIRE(r.lm.r  == Approx(1.24));
    REQUIRE(r.lm.r2 == Approx(1.24));
    REQUIRE(r.lm.r2 == Approx(1.24));
}