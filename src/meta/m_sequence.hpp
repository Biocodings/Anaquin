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
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static MSequenceStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif