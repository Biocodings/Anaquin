#ifndef V_SUBSAMPLE_HPP
#define V_SUBSAMPLE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSubsample
    {
        // How coverage is determined
        enum CoverageMethod
        {
            ArithmeticAverage,
            Maximum,
            Median,
            Percentile75,
        };

        struct Options : public AnalyzerOptions
        {
            CoverageMethod method;
        };

        struct Stats
        {
            // Coverage statistics for chrT
            Interval::Stats chrT;
            
            // Coverage statistics for human genome
            Interval::Stats hg38;

            // Raw coverage
            CoverageTool::Stats cov;
            
            // Calculated coverage for chrT
            Coverage chrTC;

            // Calculated coverage for hg38
            Coverage hg38C;

            // Fraction required to subsample in chrT
            inline double fract() const { return chrTC / hg38C; }
        };

        // Generate statistics for subsampling
        static Stats stats(const FileName &, const Options &o = Options());

        // Generate and report statistics for subsampling
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif