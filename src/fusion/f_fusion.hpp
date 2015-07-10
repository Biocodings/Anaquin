#ifndef GI_F_FUSION_HPP
#define GI_F_FUSION_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct FFusion
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : ModelStats
        {
            // Overall performance
            Performance p;

            // Overall counter
            Counter c;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif