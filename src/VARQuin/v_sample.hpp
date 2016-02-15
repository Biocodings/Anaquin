#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSample
    {
        enum CoverageMethod
        {
            ArithAverage,
            Maximum,
            Median,
            Percentile75,
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            // How coverage is calculated
            CoverageMethod method = CoverageMethod::ArithAverage;
        };

        struct Stats
        {
            // Coverage statistics for chrT
            Interval::Stats chrT;
            
            // Coverage statistics for endogenous (eg: chr21)
            Interval::Stats endo;

            // Raw coverage
            CoverageTool::Stats cov;

            // Calculated coverage for chrT
            Coverage chrTC;

            // Calculated coverage for the query (eg: chr21)
            Coverage endoC;

            /*
             * Fraction required to subsample in chrT. This works because chrT is a short
             * chromosome and almost certianly will have the highest coverage.
             */

            inline Proportion sample() const { return endoC / chrTC; }
        };

        // Generate statistics for subsampling
        static Stats stats(const FileName &, const Options &o = Options());

        // Subsample an alignment file
        static void sample(const FileName &, const FileName &, const Stats &, const Options &o = Options());
        
        // Generate and report statistics for subsampling
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
