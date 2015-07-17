#ifndef GI_F_ALIGN_HPP
#define GI_F_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FAlign
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