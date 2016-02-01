#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <data/locus.hpp>
#include <data/biology.hpp>

namespace Anaquin
{
    struct Variant
    {
        operator const Locus &() const { return l; }
        
        inline bool operator<(const Locus &x) const { return l < x; }

        // Eg: chrT
        ChromoID chrID;

        // Eg: D_1_10
        SequinID id;

        // The reference position, with the 1st base having position 1
        Locus l;
        
        // Type of the mutation
        Mutation type;
        
        Sequence ref, alt;
        
        Genotype gt;
        
        // Allelle frequency
        double af;
        
        // Allele count in genotypes
        Counts ac;
        
        // Total number of alleles in called genotypes
        Counts an;
        
        // Combined depth across samples
        unsigned dp;
        
        // Depth for reference
        unsigned dp_r;
        
        // Depth for alternative
        unsigned dp_a;
    };
}

#endif