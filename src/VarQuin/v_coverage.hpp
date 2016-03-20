#ifndef V_COVERAGE_HPP
#define V_COVERAGE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VCoverage
    {
        struct Options : public AnalyzerOptions
        {
            // The endogenous chromosome to be compared with. Optional.
            ChrID endoID;
        };

        struct Stats
        {
            // Coverage statistics for chrT
            CoverageTool::Stats chrT;

            // Coverage statistics for endogenous
            CoverageTool::Stats endo;
        };

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif