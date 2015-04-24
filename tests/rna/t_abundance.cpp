#include <catch.hpp>
#include "abundance.hpp"

using namespace Spike;

TEST_CASE("Abundance_Gene_Simulations_1")
{
    /*
    const auto r = Abundance::analyze("tests/data/rna_sims_1/assembly/genes.fpkm_tracking");

    REQUIRE(r.sb.id == "R_10_1");
    REQUIRE(r.sb.counts == 1);
    REQUIRE(r.sb.abund == 2500000);
     */
}

TEST_CASE("Abundance_Isoform_Simulations_1")
{
    /*
    Abundance::Options o;
    o.level = Abundance::Isoform;
    const auto r = Abundance::analyze("tests/data/rna_sims_1/assembly/isoforms.fpkm_tracking", o);

    REQUIRE(r.sb.id == "R_10_1_R");
    REQUIRE(r.sb.counts == 1);
    REQUIRE(r.sb.abund == 1250000);
     */
}