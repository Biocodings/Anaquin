#include <catch.hpp>
#include "rna/differential.hpp"

using namespace Spike;

TEST_CASE("Differential_Gene_RNA_Sims_2")
{
    const auto r = Differential::analyze("tests/data/rna_sims_1/diffs/gene_exp.diff");
}

TEST_CASE("Differential_Isoforms_RNA_Sims_2")
{
    Differential::Options o;
    o.mode = Differential::DiffIsoform;
    const auto r = Differential::analyze("tests/data/rna_sims_1/diffs/isoform_exp.diff", o);
}