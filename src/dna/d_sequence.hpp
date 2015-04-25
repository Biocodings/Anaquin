#ifndef GI_D_SEQUENCE_HPP
#define GI_D_SEQUENCE_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct DSequenceStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct DSequence
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static DSequenceStats analyze(const std::string &file, const DSequence &options = DSequence());
    };
}

#endif