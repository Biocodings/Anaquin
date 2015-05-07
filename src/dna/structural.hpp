#ifndef GI_STRUCTURAL_HPP
#define GI_STRUCTURAL_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct StructuralStats
    {
        
        
        
    };

    struct Structural
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static StructuralStats analyze(const std::string &file, const Options &options = Options());
    };    
}

#endif