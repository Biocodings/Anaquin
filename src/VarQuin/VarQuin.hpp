#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/data.hpp"
#include "data/variant.hpp"

namespace Anaquin
{
    /*
     * Information specific to sequin variants.
     */
    
    struct SeqVariant
    {
        enum class Group
        {
            NA12878,
            VeryLowGC,
            LowGC,
            HighGC,
            VeryHighGC,
            ShortDinRep,  // Dinucleotide repeats
            LongDinRep,   // Dinucleotide repeats
            ShortHompo,
            LongHompo,
            ShortQuadRep, // Quad-nucleotide repeats
            LongQuadRep,  // Quad-nucleotide repeats
            ShortTrinRep, // Trinucleotide repeats
            LongTrinRep,  // Trinucleotide repeats
            Cosmic,
        } group;
        
        // Germline? Somatic?
        Mutation mut;
        
        // Homozygous?
        Zygosity zyg;
        
        // Copy number
        unsigned copy = 1;
        
        // Valid only for group == Cosmis
        std::string info;
    };
    
    // Eg: chrev1, chrev10 etc...
    inline bool isReverseGenome(const ChrID &cID)
    {
        A_ASSERT(!cID.empty());
        return cID.find("rev") != std::string::npos;
    }
}

#endif
