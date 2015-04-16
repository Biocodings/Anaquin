#ifndef GI_ALIGNER_HPP
#define GI_ALIGNER_HPP

#include "analyzer.hpp"
#include "confusion.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct AlignerStats : public AnalyzerStats
    {
        // The lowest detectable abundance in the experiment
        Sensitivity sens;

        Confusion m;
    };

    struct Aligner
    {
        enum Mode
        {
            AlignBase,
            AlignExon,
            AlignSplice,
        };

        struct Options : public AnalyzerOptions<Aligner::Mode>
        {
            // Empty Implementation
        };

        static AlignerStats analyze(const std::string &file, const Aligner::Options &options = Aligner::Options());
    };
}

#endif