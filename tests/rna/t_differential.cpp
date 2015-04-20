#include <catch.hpp>
#include "rna/differential.hpp"

using namespace Spike;

TEST_CASE("Differential_Genes_Simulations_1")
{
    const auto r = Differential::analyze("tests/data/rna_sims_1/diffs/gene_exp.diff");
}

TEST_CASE("Differential_Isoforms_Simulations_1")
{
    Differential::Options o;
    o.level = Differential::Isoform;
    const auto r = Differential::analyze("tests/data/rna_sims_1/diffs/isoform_exp.diff", o);
}