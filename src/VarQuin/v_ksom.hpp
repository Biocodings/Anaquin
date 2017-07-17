#ifndef V_KSOMATIC_HPP
#define V_KSOMATIC_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VarKSomatic
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
