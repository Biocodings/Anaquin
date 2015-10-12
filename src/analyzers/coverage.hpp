#ifndef COVERAGE_ANALYZER_HPP
#define COVERAGE_ANALYZER_HPP

#include "stats/analyzer.hpp"
#include "tools/coverage.hpp"

namespace Anaquin
{
    class CoverageAnalyzer : public Analyzer
    {
        public:
        
            typedef CoverageTool::Stats Stats;
            typedef AnalyzerOptions Options;

            static CoverageTool::Stats stats(const FileName &, const Options &o = Options());
            static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif