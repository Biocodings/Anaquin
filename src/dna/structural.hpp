#ifndef GI_STRUCTURAL_HPP
#define GI_STRUCTURAL_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct StructuralStats
    {
        // Overall performance
        Performance p;

        // Performance relative to alleles
        Performance p_al;
        
        // Performance relative to the position
        Performance p_l;

        // Performance relative to genotype
        Performance p_gt;

        // Performance relative to allele frequency
        Performance p_af;
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