#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <iostream>
#include <data/locus.hpp>
#include <data/biology.hpp>

namespace Anaquin
{
    struct Variant
    {
        operator const Locus &() const { return l; }
        
        inline bool operator<(const Locus &x) const { return l < x; }

        inline Mutation type() const
        {
            if (alt[0] == '-')
            {
                return Deletion;
            }
            else if (alt[0] == '+')
            {
                return Insertion;
            }
            else if (ref.size() == alt.size())
            {
                return SNP;
            }
            else if (ref.size() > alt.size())
            {
                return Deletion;
            }
            else
            {
                return Insertion;
            }
        }

        inline Base diff() const
        {
            return labs((Base)ref.size() - (Base)alt.size());
        }

        // Eg: chrT
        ChrID chrID;

        // Eg: D_1_10
        SequinID id;

        // The reference position, with the 1st base having position 1
        Locus l;
        
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
    
    struct CalledVariant : public Variant
    {
        inline Proportion alleleFreq() const
        {
            if (readR && readV)
            {
                return static_cast<Proportion>(readV) / (readR + readV);
            }

            return static_cast<Proportion>(dp_a) / dp_r;
        }
    
        Probability p;

        // Number of reads for the reference
        Counts readR = 0;
        
        // Number of reads for the variant
        Counts readV = 0;
    };
}

#endif