#ifndef GI_C_DIFFS_HPP
#define GI_C_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LDiffs
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