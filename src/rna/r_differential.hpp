#ifndef GI_R_DIFFERENTIAL_HPP
#define GI_R_DIFFERENTIAL_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RDifferentialStats : public AnalyzerStats
    {
        Confusion m;
        Sensitivity s;

        // Correlation for the samples
        double r;

        // Adjusted R2 for the linear model
        double r2;

        // Coefficient for the linear model
        double slope;
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