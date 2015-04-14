#include <catch.hpp>
#include "abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_RNA_Simulation_2")
{
    const auto r = Abundance::analyze("tests/data/rna_sims_2/assembled/transcripts.gtf");

}