#include <catch.hpp>
#include "rna/r_abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_Gene_Simulations_1")
{
    const auto r = RAbundance::analyze("tests/data/rna_sims/genes.fpkm_tracking");

    REQUIRE(r.sb.id == "R_10_1");
    REQUIRE(r.sb.counts == 1);
    REQUIRE(r.sb.abund == 2500000);
}

TEST_CASE("Abundance_Isoform_Simulations_1")
{
    //RAbundance::Options o;
    //o.level = RAbundance::Isoform;
    //const auto r = RAbundance::analyze("tests/data/rna_sims/isoforms.fpkm_tracking", o);

    //REQUIRE(r.sb.id == "R_10_1_R");
    //REQUIRE(r.sb.counts == 1);
    //REQUIRE(r.sb.abund == 1250000);
}