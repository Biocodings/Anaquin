#ifndef GI_ALIGNER_HPP
#define GI_ALIGNER_HPP

#include "analyzer.hpp"
#include "classify.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct AlignerStats : public AnalyzerStats
    {
        // The lowest detectable abundance in the experiment
        Sensitivity sens;

        SS::Confusion m;
    };

    struct Aligner
    {
        enum Mode
        {
            Base,
            Exon,
            Splice,
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