#ifndef GI_STRUCTURAL_HPP
#define GI_STRUCTURAL_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct StructuralStats
    {
        // Overall performance
        Confusion m;

        // Overall sensitivity
        Sensitivity s;

        // Performance for the position
        Confusion ml;

        // Sensitivity for the position
        Sensitivity sl;
    };

    struct Structural
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static StructuralStats analyze(const std::string &file, const Options &options = Options());
    };    
}

#endif