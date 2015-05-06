#ifndef GI_R_DIFFERENTIAL_HPP
#define GI_R_DIFFERENTIAL_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RDifferentialStats : public AnalyzerStats
    {
        Sensitivity s;
        LinearModel lm;
    };

    struct RDifferential : public RAnalyzer
    {
        struct Options : public DoubleMixtureOptions
        {
            RNALevel level;
        };

        static RDifferentialStats analyze(const std::string &f, const Options &options = Options());
    };
}

#endif