#ifndef GI_D_SEQUENCE_HPP
#define GI_D_SEQUENCE_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct DSequenceStats
    {
        // Empty Implementation
    };

    struct DSequence
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static DSequenceStats analyze(const std::string &file, const DSequence &options = DSequence());
    };
}

#endif