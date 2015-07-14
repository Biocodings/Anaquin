#ifndef GI_R_ABUND_HPP
#define GI_R_ABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RAbundStats : public ModelStats
    {
        // Empty Implementation
    };

    struct RAbund : public RAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            Options() {}
            RNALevel level = Isoform;
        };

        static RAbundStats analyze(const std::string &, const Options &options = Options());
    };
}

#endif