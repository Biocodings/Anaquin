#ifndef GI_L_NORM_HPP
#define GI_L_NORM_HPP

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

        static Stats report(const std::string &, const std::string &, const Options &options = Options());
    };
}

#endif