#include <catch.hpp>
#include "abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_Gene_RNA_Sims_2")
{
    const auto r = Abundance::analyze("tests/data/rna_sims_2/assembly/genes.fpkm_tracking",
                      Abundance::AbundanceOptions(Abundance::AbdunanceMode::AbdunanceGene));
}

TEST_CASE("Abundance_Isoform_RNA_Sims_2")
{
    const auto r = Abundance::analyze("tests/data/rna_sims_2/assembly/isoforms.fpkm_tracking",
                      Abundance::AbundanceOptions(Abundance::AbdunanceMode::AbdunanceIsoform));
}