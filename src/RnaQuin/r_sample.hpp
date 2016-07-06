#ifndef R_SAMPLE_HPP
#define R_SAMPLE_HPP

#include "tools/sample.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RSample
    {
        enum class Method
        {
            Prop
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            // Defined if the method is proportion
            Proportion p;

            Method meth = Method::Prop;
        };

        struct Stats : public MappingStats
        {
            typedef Sampler::SGReads SGReads;
            
            // Reads before and after subsampling
            SGReads before, after;

            // Proportion required for subsampling
            Proportion p = NAN;
        };

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif