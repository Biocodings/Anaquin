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
        enum DifferentialMode
        {
            DiffGene,
            DiffIsoform
        };

        struct DifferentialOptions : public AnalyzerOptions
        {
            DifferentialOptions(DifferentialMode mode) : mode(mode) {}
            
            // Whether it's done at the gene or isoform level
            DifferentialMode mode;
        };

        static DifferentialStats analyze(const std::string &f, const DifferentialOptions &options);
    };
}

#endif