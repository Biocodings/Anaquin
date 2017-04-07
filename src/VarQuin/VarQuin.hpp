#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/data.hpp"
#include "data/variant.hpp"

namespace Anaquin
{
    struct VariantStats
    {
        // Number of SNPs detected
        Counts n_snp;

        // Number of indels detected
        Counts n_ind;
    };
    
    struct VariantMatch
    {
        // The called variant
        Variant query;

        // Matched by position?
        const Variant *match = nullptr;
        
        // Matched by variant allele? Only if position is matched.
        bool alt;
        
        // Matched by reference allele? Only if position is matched.
        bool ref;
    };

    // Eg: chrev1, chrev10 etc...
    inline bool isReverseGenome(const ChrID &cID)
    {
        A_ASSERT(!cID.empty());
        return cID.find("rev") != std::string::npos;
    }
}

#endif
