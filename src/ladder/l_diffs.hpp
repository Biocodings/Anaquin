#ifndef L_DIFFS_HPP
#define L_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LDiffs
    {
        typedef AnalyzerOptions Options;
        typedef LinearStats Stats;

        static Stats report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif