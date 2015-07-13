#ifndef GI_R_ABUNDANCE_HPP
#define GI_R_ABUNDANCE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RAbundanceStats : public ModelStats
    {
        // Empty Implementation
    };

    struct RAbundance : public RAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            Options() {}
            RNALevel level = Isoform;
        };

        static RAbundanceStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif