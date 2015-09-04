#ifndef GI_T_ASSEMBLY_HPP
#define GI_T_ASSEMBLY_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TAssembly : Analyzer
    {
        struct Options : FuzzyOptions
        {
            // Path for the reference and query GTF
            std::string ref, query;
        };
        
        struct Stats
        {
            BaseHist   hb = Analyzer::baseHist();
            SequinHist he = Analyzer::seqHist();
            SequinHist hi = Analyzer::seqHist();
            SequinHist ht = Analyzer::seqHist();

            // Overall performance
            Performance p;
            
            // Performance at the base level
            Performance pb;
            
            // Performance at the exon level
            Performance pe;
            
            // Performance at the transcript level
            Performance pt;
            
            // Performance at the intron level
            Performance pi;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif