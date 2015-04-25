#ifndef GI_R_SEQUENCE_HPP
#define GI_R_SEQUENCE_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RSequenceStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct RSequence : RAnalyzer
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static RSequenceStats analyze(const std::string &file, const RSequence &options = RSequence());
    };
}

#endif