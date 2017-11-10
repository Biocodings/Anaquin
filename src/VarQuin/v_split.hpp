#ifndef V_SPLIT_HPP
#define V_SPLIT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSplit
    {
        typedef AnalyzerOptions Options;

        struct Stats
        {
            Counts nNMap  = 0;
            Counts nAmbig = 0;
            
            // Number of sample-derived alignment reads
            Counts nSample = 0;
            
            // Number of germline alignment reads
            Counts nGerm = 0;
            
            // Number of somatic alignment reads
            Counts nSom = 0;
            
            // Number of strucutral variant alignment reads
            Counts nStru = 0;
            
            // Number of ladder alignment reads
            Counts nLad = 0;
        };

        static Stats analyze(const FileName &, const Options &);
        static void report(const FileName &, const Options &);
    };
}

#endif
