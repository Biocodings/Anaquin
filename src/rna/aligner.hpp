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
            LevelBase,
            LevelExon,
            LevelSplice,
        };

        struct Options : public AnalyzerOptions<Aligner::Mode>
        {
            // Empty Implementation
        };

        inline static std::string name() { return "align"; }

        static AlignerStats analyze(const std::string &file, const Aligner::Options &options = Aligner::Options());
    };
}

#endif