#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"
#include "tools/subsample.hpp"

namespace Anaquin
{
    struct VSample
    {
        struct Options : public AnalyzerOptions
        {
            Options() {}

            // How coverage is calculated
            Subsampler::CoverageMethod method = Subsampler::CoverageMethod::ArithAverage;
        };

        typedef Subsampler::Stats Stats;

        static Stats stats(const FileName &file, const Options &o = Options())
        {
            return Subsampler::stats(file, o);
        }

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
