#ifndef GI_D_SEQUENCE_HPP
#define GI_D_SEQUENCE_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct DSequenceStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct DSequence
    {
        enum Mode
        {
            Base,
        };

        struct Options : public AnalyzerOptions<DSequence::Mode>
        {
            // Empty Implementation
        };

        static DSequenceStats analyze(const std::string &file, const DSequence::Options &options = DSequence::Options());
    };
}

#endif