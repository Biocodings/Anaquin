#ifndef GI_RALIGNER_HPP
#define GI_RALIGNER_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct RAlignerStats : public AnalyzerStats
    {
        Sensitivity se;
    };

    struct RAligner
    {
        enum Mode
        {
            Base,
            Exon,
            Splice,
        };

        struct Options : public AnalyzerOptions<RAligner::Mode>
        {
            // Empty Implementation
        };

        static RAlignerStats analyze(const std::string &file, const RAligner::Options &options = RAligner::Options());
    };
}

#endif