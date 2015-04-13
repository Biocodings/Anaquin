#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include "types.hpp"
#include "parsers/parser.hpp"

namespace Spike
{
    struct AbundanceStats
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
        struct AbundanceOptions : public ParserOptions
        {
            // Empty Implementation
        };
        
        static AbundanceStats analyze(const std::string &file, const AbundanceOptions &options = AbundanceOptions());
    };    
}

#endif