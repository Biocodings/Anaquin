#ifndef GI_R_ABUNDANCE_HPP
#define GI_R_ABUNDANCE_HPP

#include "r_analyzer.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct RAbundanceStats : public AnalyzerStats
    {
        // Correlation for the samples
        double r;
        
        // Adjusted R2 for the linear model
        double r2;
        
        // Coefficient for the linear model
        double slope;
    };

    struct RAbundance : public RAnalyzer
    {
        enum Level
        {
            Gene,
            Isoform,
        };

        struct Options : public AnalyzerOptions<RAbundance::Level>
        {
            // Empty Implementation
        };

        static RAbundanceStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif