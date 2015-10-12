#ifndef V_COVERAGE_HPP
#define V_COVERAGE_HPP

#include "analyzers/coverage.hpp"

namespace Anaquin
{
    struct VCoverage
    {
        typedef AnalyzerOptions Options;
        typedef CoverageAnalyzer::Stats Stats;
        
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif