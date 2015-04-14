#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include "analyzer.hpp"
#include "confusion.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct AbundanceStats : public AnalyzerStats
    {
        Confusion m_base_;

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
        enum AbdunanceMode
        {
            AbdunanceGene,
            AbdunanceIsoform
        };
        
        struct AbundanceOptions : public AnalyzerOptions
        {
            AbundanceOptions(AbdunanceMode mode) : mode(mode) {}

            AbdunanceMode mode;
        };

        static AbundanceStats analyze(const std::string &file, const AbundanceOptions &options);
    };    
}

#endif