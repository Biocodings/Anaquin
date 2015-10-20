#ifndef L_COVERAGE_HPP
#define L_COVERAGE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LCoverage
    {
        typedef AnalyzerOptions Options;        
        typedef CoverageTool::Stats Stats;

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif