//#include <catch.hpp>
//#include "rna/r_differential.hpp"
//
//using namespace Spike;
//
//TEST_CASE("Differential_Simulations_Isoforms")
//{
//    /*
//     * The experiment for mixutre B is identical but it has 20 times coverage.
//     */
//    
//    for (int i = 0; i < 5000; i++)
//    {
//    
//    RDifferential::Options o;
//    o.level = RNALevel::Isoform;
//    const auto r = RDifferential::analyze("tests/data/rna_sims/isoform_exp.diff", o);
//    }
//    
//
//}
//
//TEST_CASE("Differential_Simulations_Genes")
//{
//    /*
//     * The experiment for mixutre B is identical but it has 20 times coverage.
//     */
//    
//    RDifferential::Options o;
//    o.level = RNALevel::Gene;
//    const auto r = RDifferential::analyze("tests/data/rna_sims/gene_exp.diff", o);
//
//    
//
//}