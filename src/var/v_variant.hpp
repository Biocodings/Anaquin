#ifndef GI_D_VARIANT_HPP
#define GI_D_VARIANT_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct DVariant
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public ModelStats
        {
            // Overall performance
            Confusion m;

            // The proportion of variations with alignment coverage
            double covered;

            // Measure of variant detection independent to sequencing depth or coverage
            double efficiency;

            SequinCounter c = DAnalyzer::counterSequins();
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif