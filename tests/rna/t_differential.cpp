#include <catch.hpp>
#include "rna/r_differential.hpp"

using namespace Spike;

TEST_CASE("Differential_Genes_Simulations_1")
{
    const auto r = RDifferential::analyze("tests/data/rna_sims/gene_exp.diff");
}

TEST_CASE("Differential_Isoforms_Simulations_1")
{
    //RDifferential::Options o;
    //o.level = RDifferential::Isoform;
    //const auto r = RDifferential::analyze("tests/data/rna_sims/isoform_exp.diff", o);
}