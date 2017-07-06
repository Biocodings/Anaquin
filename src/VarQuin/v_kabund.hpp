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
             * Germline mutation ladder
             */
            
            std::map<SequinID, Measured> germR, germV;
            
            /*
             * Cancer mutation ladder
             */
            
            std::map<SequinID, Measured> canR, canV;

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
