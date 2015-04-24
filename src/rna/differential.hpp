#ifndef GI_DIFFERENTIAL_HPP
#define GI_DIFFERENTIAL_HPP

#include "classify.hpp"
#include "r_analyzer.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct DifferentialStats : public AnalyzerStats
    {
        // Correlation for the samples
        double r;

        // Adjusted R2 for the linear model
        double r2;

        // Coefficient for the linear model
        double slope;
    };

    struct Differential : public RAnalyzer
    {
        enum Level
        {
            Gene,
            Isoform,
        };

        struct Options : public AnalyzerOptions<Differential::Level>
        {
            // Empty Implementation
        };

        static DifferentialStats analyze(const std::string &f, const Options &options = Options());
    };
}

#endif