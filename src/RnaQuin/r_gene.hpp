#ifndef R_GENE_HPP
#define R_GENE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RGene
    {
        typedef AnalyzerOptions Options;
        
        struct Data
        {
            GeneID gID;
            
            // Length of the gene
            Base len;
            
            // Concentration at the gene level
            Concent con;
        };

        typedef std::map<Mixture, std::map<GeneID, Data>> Stats;

        static Stats stats(const Options &o);
        static void report(const Options &o = Options());
    };
}

#endif