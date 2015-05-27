#ifndef GI_D_VARIANT_HPP
#define GI_D_VARIANT_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct VariantStats
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

    struct DVariant
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static VariantStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif