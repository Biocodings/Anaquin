#ifndef GI_R_DIFFS_HPP
#define GI_R_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct RDiffs : public RAnalyzer
    {
        struct Options : public DoubleMixtureOptions
        {
            RNALevel level;
        };

        struct Stats : public ModelStats
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &f, const Options &options = Options());
    };
}

#endif