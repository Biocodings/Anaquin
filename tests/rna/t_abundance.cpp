#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_Gene_Simulations_1")
{
    RAbundance::Options o;
    o.level = RNALevel::Gene;
  //const auto r = RAbundance::analyze("tests/data/rna_sims/genes.fpkm_tracking", o);
}

TEST_CASE("Abundance_Isoform_Simulations_1")
{
    RAbundance::Options o;
    o.level = RNALevel::Isoform;
   // const auto r = RAbundance::analyze("tests/data/rna_sims/isoforms.fpkm_tracking", o);
}