#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include "data/locus.hpp"

namespace Anaquin
{
    struct Alignment
    {
        operator const Locus &() const { return l; }

        // Eg: B7_591:6:155:12:674
        ReadID name;

        // Primary alignment
        ChrID cID;

        // Location of the alignment
        Locus l;
        
        // If this field is false, no assumption can be made to other fields
        bool mapped;
        
        // Only valid if the alignment is spliced
        Base skipped;

        // Mapping quality
        int mapq;

        // Bitwise FLAG
        int flag;

        /*
         * SAM flag fields
         */
        
        bool isPaired;
        bool isAllAligned;
        bool isAligned;
        bool isMateAligned;
        bool isForward;
        bool isMateReverse;
        bool isFirstPair;
        bool isSecondPair;
        bool isDuplicate;
        bool isPrimary;
        bool isSupplement;        
        bool isPassed;
        
        // Secondary alignment? Typically used for alternative mappings when multiple mappings are presented
        bool isSecondary;
        
        /*
         * Optional fields
         */
        
        // Signed observed template length
        Base tlen;
        
        // Cigar string
        std::string cigar;

        // Reference sequence name of the primary alignment
        ChrID rnext;
        
        // Position of the primary alignment of the NEXT read in the template
        std::string pnext;

        // Segment sequence
        std::string seq;
        
        // ASCII of base QUALity plus 33
        std::string qual;        
    };
}

#endif
