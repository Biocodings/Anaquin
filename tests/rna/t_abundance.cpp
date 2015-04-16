#include <catch.hpp>
#include "abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_Gene_RNA_Sims_2")
{
    const auto r = Abundance::analyze("tests/data/rna_sims_1/assembly/genes.fpkm_tracking");
}

TEST_CASE("Abundance_Isoform_RNA_Sims_2")
{
    Abundance::Options o;
    o.mode = Abundance::AbdunanceIsoform;
    const auto r = Abundance::analyze("tests/data/rna_sims_1/assembly/isoforms.fpkm_tracking", o);
}