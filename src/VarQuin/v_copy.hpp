#ifndef V_COPY_HPP
#define V_COPY_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VCopy
    {
        typedef unsigned CopyNumber;
        
        struct Options : public AnalyzerOptions
        {
            Options() : copy(2) {}
            
            CopyNumber copy;
        };
        
        struct Stats
        {
        };

        static Stats analyze(const FileName &, const Options &o);
        static void report  (const FileName &, const Options &o = Options());
    };
}

#endif
