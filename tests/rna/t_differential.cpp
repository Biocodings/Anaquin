#include <catch.hpp>
#include "rna/differential.hpp"

using namespace Spike;

TEST_CASE("Differential_Gene_RNA_Sims_2")
{
    const auto r = Differential::analyze("tests/data/rna_sims_2/diffs/gene_exp.diff", Differential::DifferentialOptions(Differential::DiffGene));
}

TEST_CASE("Differential_Isoforms_RNA_Sims_2")
{
    const auto r = Differential::analyze("tests/data/rna_sims_2/diffs/isoform_exp.diff", Differential::DifferentialOptions(Differential::DiffIsoform));
}