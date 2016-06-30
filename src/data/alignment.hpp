#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include "data/locus.hpp"

namespace Anaquin
{
    struct Alignment
    {
        operator const Locus &() const { return l; }

        ChrID cID;

        // Location of the alignment
        Locus l;
        
        // If this field is false, no assumption can be made to other fields
        bool mapped;
        
        // Only valid if the alignment is spliced
        Base skipped;

        /*
         * Optional fields
         */
        
        // Eg: B7_591:6:155:12:674
        std::string name;

        // Signed observed template length
        Base tlen;
        
        // Cigar string
        std::string cigar;

        // Bitwise FLAG
        int flag;

        // Reference sequence name of the primary alignment
        std::string rnext;
        
        // Position of the primary alignment of the NEXT read in the template
        std::string pnext;

        // Mapping quality
        int mapq;
        
        // Segment sequence
        std::string seq;
        
        // ASCII of base QUALity plus 33
        std::string qual;        
    };
}

#endif