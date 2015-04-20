#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include "analyzer.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct AbundanceStats : public AnalyzerStats
    {
        Sensitivity s;
        
        // Correlation for the samples
        double r;
        
        // Adjusted R2 for the linear model
        double r2;
        
        // Coefficient for the linear model
        double slope;
    };
    
    struct Abundance
    {
        struct Options : public AnalyzerOptions<RNALevel>
        {
            // Empty Implementation
        };

        static AbundanceStats analyze(const std::string &file, const Abundance::Options &options = Abundance::Options());
    };
}

#endif