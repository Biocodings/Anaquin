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

        struct Stats
        {

        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif