#include <catch.hpp>
#include "rna/r_differential.hpp"

using namespace Spike;

TEST_CASE("Differential_Simulations_Isoforms")
{
    RDifferential::Options o;
    o.level = RNALevel::Isoform;
 
    
    const auto r = RDifferential::analyze("tests/data/rna_sims/isoform_exp.diff", o);


}

TEST_CASE("Differential_Simulations_Genes")
{
    RDifferential::Options o;
    o.level = RNALevel::Gene;
    
    
    const auto r = RDifferential::analyze("tests/data/rna_sims/gene_exp.diff", o);


}