#ifndef GI_DNA_VARIATION_HPP
#define GI_DNA_VARIATION_HPP

#include "types.hpp"
#include "analyzer.hpp"

namespace Spike
{
    struct VariationStats : public AnalyzerStats
    {
        Percentage covered;
        Percentage efficiency;
    };

    struct DNAVariation
    {
        static VariationStats analyze(const std::string &file);
    };    
}

#endif