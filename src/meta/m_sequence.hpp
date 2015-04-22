#ifndef GI_M_SEQUENCE_HPP
#define GI_M_SEQUENCE_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct MSequenceStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct MSequence
    {
        enum Mode
        {
            Base,
        };

        struct Options : public AnalyzerOptions<MSequence::Mode>
        {
            // Empty Implementation
        };

        static MSequenceStats analyze(const std::string &file, const MSequence &options = MSequence());
    };
}

#endif