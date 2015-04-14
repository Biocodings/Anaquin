#ifndef GI_DIFFERENTIAL_HPP
#define GI_DIFFERENTIAL_HPP

#include "analyzer.hpp"
#include "confusion.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct DifferentialStats : public AnalyzerStats
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

    struct Differential
    {
        struct DifferentialOptions : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static DifferentialStats analyze(const std::string &s1,
                                         const std::string &s2,
                                         const DifferentialOptions &options);
    };
}

#endif