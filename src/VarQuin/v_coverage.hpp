#ifndef V_COVERAGE_HPP
#define V_COVERAGE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VCoverage
    {
        typedef AnalyzerOptions Options;
        typedef CoverageTool::Stats Stats;

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif