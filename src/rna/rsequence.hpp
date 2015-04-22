#ifndef GI_R_SEQUENCE_HPP
#define GI_R_SEQUENCE_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct RSequenceStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct RSequence
    {
        enum Level
        {
            Base
        };

        struct Options : public AnalyzerOptions<RSequence::Level>
        {
            // Empty Implementation
        };

        static RSequenceStats analyze(const std::string &file, const RSequence::Options &options = RSequence::Options());
    };
}

#endif