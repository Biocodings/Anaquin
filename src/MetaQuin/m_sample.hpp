#ifndef M_SAMPLE_HPP
#define M_SAMPLE_HPP

#include "tools/sample.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MSample
    {
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Subsampler::CoverageMethod method = Subsampler::CoverageMethod::ArithAverage;
        };

        typedef Subsampler::Stats Stats;

        static Stats stats(const FileName &file, const Options &o = Options())
        {
            return Subsampler::stats(file, o);
        }

        static void report(const FileName &file, const Options &o = Options())
        {
            Subsampler::report(file, o);
        }
    };
}

#endif