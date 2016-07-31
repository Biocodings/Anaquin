#ifndef R_SAMPLE_HPP
#define R_SAMPLE_HPP

#include "tools/sample.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RSample
    {
        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            bool toConsole = true;
            
            // Fraction required for the spike-in
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

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif