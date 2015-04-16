#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include "analyzer.hpp"
#include "confusion.hpp"
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
        enum Mode
        {
            AbdunanceGene,
            AbdunanceIsoform
        };

        struct Options : public AnalyzerOptions<Abundance::Mode>
        {
            // Empty Implementation
        };

        inline static std::string name() { return "abundance"; }

        static AbundanceStats analyze(const std::string &file, const Abundance::Options &options = Abundance::Options());
    };
}

#endif