#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include "analyzer.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct AbundanceStats : public AnalyzerStats
    {
        // Correlation for the samples
        double r;
        
        // Adjusted R2 for the linear model
        double r2;
        
        // Coefficient for the linear model
        double slope;
    };
    
    struct Abundance
    {
        enum Level
        {
            Gene,
            Isoform,
        };

        struct Options : public AnalyzerOptions<Abundance::Level>
        {
            // Empty Implementation
        };

        static AbundanceStats analyze(const std::string &file, const Abundance::Options &options = Abundance::Options());
    };
}

#endif