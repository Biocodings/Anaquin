#ifndef GI_F_DISCOVER_HPP
#define GI_F_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FDiscover
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats
        {
            // Overall performance
            Confusion m;

            // Distribution of the sequins
            SequinHist h;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif