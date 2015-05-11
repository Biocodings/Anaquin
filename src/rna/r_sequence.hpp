#ifndef GI_R_SEQUENCE_HPP
#define GI_R_SEQUENCE_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct RSequenceStats
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