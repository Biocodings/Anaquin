#ifndef V_KABUND_HPP
#define V_KABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKAbund
    {
        struct Stats
        {
            // Read counts for reference and variants (allele frequency ladder)
            std::map<SequinID, Measured> afR, afV;
            
            // Statistics for conjoint sequins
            SequinStats con;            
        };
        
        typedef AnalyzerOptions Options;
        
        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
