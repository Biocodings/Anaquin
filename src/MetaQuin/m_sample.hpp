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
            
            // Fraction required for sampling
            Proportion p = NAN;
        };

        struct Stats : public MappingStats
        {
            typedef Sampler::SGReads SGReads;
            
            // Reads before and after subsampling
            SGReads before, after;

            // Normalization factor
            Proportion norm;
        };

        static Stats analyze(const FileName &, const Options &);
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
