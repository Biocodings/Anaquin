#ifndef V_KABUND_HPP
#define V_KABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKAbund
    {
        struct Stats
        {
            /*
             * Copy number ladder
             */
            
            /*
             * Allele frequency ladder
             */
            
            std::map<SequinID, Measured> afR, afV;

            /*
             * Conjoint ladder
             */
            
            std::map<SequinID, Measured> con;
        };
        
        typedef AnalyzerOptions Options;
        
        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
