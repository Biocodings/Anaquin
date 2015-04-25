#ifndef GI_R_DIFFERENTIAL_HPP
#define GI_R_DIFFERENTIAL_HPP

#include "classify.hpp"
#include "r_analyzer.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct RDifferentialStats : public AnalyzerStats
    {
        // Correlation for the samples
        double r;

        // Adjusted R2 for the linear model
        double r2;

        // Coefficient for the linear model
        double slope;
    };

    struct RDifferential : public RAnalyzer
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static RDifferentialStats analyze(const std::string &f, const Options &options = Options());
    };
}

#endif