#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include "analyzer.hpp"
#include "sensitivity.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    struct AbundanceStats : public AnalyzerStats
    {
        ConfusionMatrix m_base;

        Sensitivity s_base;
        
        // Correlation for the samples
        double r;
        
        // Adjusted R2 for the linear model
        double r2;
        
        // Coefficient for the linear model
        double slope;
    };
    
    struct Abundance
    {
        struct AbundanceOptions : public AnalyzerOptions
        {
            // Empty Implementation
        };
        
        static AbundanceStats analyze(const std::string &file, const AbundanceOptions &options = AbundanceOptions());
    };    
}

#endif