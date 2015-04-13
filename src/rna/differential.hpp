#ifndef GI_DIFFERENTIAL_HPP
#define GI_DIFFERENTIAL_HPP

#include "types.hpp"
#include "parsers/parser.hpp"

namespace Spike
{
    struct DifferentialStats
    {
        // Correlation for the samples
        double r;
        
        // Adjusted R2 for the linear model
        double r2;
        
        // Coefficient for the linear model
        double slope;
    };
    
    struct Differential
    {
        struct DifferentialOptions : public ParserOptions
        {
            // Empty Implementation
        };

        static DifferentialStats analyze(const std::string &file, const DifferentialOptions &options = DifferentialOptions());
    };    
}

#endif