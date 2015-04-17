#include <catch.hpp>
#include "abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_Gene_1")
{
    const auto r = Abundance::analyze("tests/data/rna_sims_1/assembly/genes.fpkm_tracking");
}

TEST_CASE("Abundance_Isoform_1")
{
    Abundance::Options o;
    o.level = LevelIsoform;
    const auto r = Abundance::analyze("tests/data/rna_sims_1/assembly/isoforms.fpkm_tracking", o);
}