#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include "data/locus.hpp"

namespace Anaquin
{
    struct Alignment
    {
        operator const Locus &() const { return l; }

        ChrID cID;

        // Eg: B7_591:6:155:12:674
        std::string name;

        // Bitwise FLAG
        int flag;
        
        // Mapping quality
        int mapq;
        
        // Cigar string
        std::string cigar;
        
        // Reference sequence name of the primary alignment
        std::string rnext;
        
        // Position of the primary alignment of the NEXT read in the template
        Base pnext;
        
        // Signed observed template length
        Base tlen;
        
        // Segment sequence
        std::string seq;
        
        // ASCII of base QUALity plus 33
        std::string qual;
        
        // Location of the alignment
        Locus l;

        unsigned i;
        
        // If this field is false, no assumption can be made to other fields
        bool mapped;

        // Whether this is a spliced read
        bool spliced;
        
        // Only valid if the alignment is spliced
        Base skipped;
    };
}

#endif