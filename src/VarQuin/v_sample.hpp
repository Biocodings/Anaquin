#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "tools/sample.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSample
    {
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Subsampler::CoverageMethod method = Subsampler::CoverageMethod::ArithAverage;
        };

        typedef Subsampler::Stats Stats;

        static Stats stats(const FileName &file, const Options &o = Options());
        static void report(const FileName &file, const Options &o = Options());
    };
}

#endif