#ifndef GI_L_NORM_HPP
#define GI_L_NORM_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LNorm
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public LinearStats
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &, const std::string &, const Options &options = Options());
    };
}

#endif