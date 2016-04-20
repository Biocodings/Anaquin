#ifndef M_COVERAGE_HPP
#define M_COVERAGE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MCoverage
    {
        typedef AnalyzerOptions Options;
        typedef CoverageTool::Stats Stats;
        
        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif