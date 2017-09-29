#ifndef V_KMER_HPP
#define V_KMER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKmer
    {
        struct Stats : public SequinStats
        {
            // Ladder for reference and variant sequins
            std::map<SequinID, Measured> r, v;
        };
        
        typedef AnalyzerOptions Options;
        
        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
