#ifndef L_NORM_HPP
#define L_NORM_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LNorm
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public LinearStats
        {
            // Empty Implementation
        };

        static Stats report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif