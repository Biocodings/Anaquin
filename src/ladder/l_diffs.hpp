#ifndef GI_C_DIFFS_HPP
#define GI_C_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LDiffs
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public ModelStats
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &, const std::string &, const Options &options = Options());
    };
}

#endif