#ifndef GI_R_ABUNDANCE_HPP
#define GI_R_ABUNDANCE_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct RAbundanceStats : public CorrelationStats
    {
        // Empty Implementation
    };

    struct RAbundance : public RAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            RNALevel level = Isoform;
        };

        static RAbundanceStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif