#ifndef F_COVERAGE_HPP
#define F_COVERAGE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FCoverage
    {
        typedef AnalyzerOptions Options;        
        typedef CoverageTool::Stats Stats;

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif