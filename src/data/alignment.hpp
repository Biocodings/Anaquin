#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include "data/locus.hpp"

namespace Anaquin
{
    struct Alignment
    {
        typedef std::string AlignmentID;

        operator const Locus &() const { return l; }

        AlignmentID id;

        // Eg: B7_591:6:155:12:674
        std::string name;
        
        unsigned i;
        
        // If this field is false, no assumption can be made to other fields
        bool mapped;

        // Whether this is a spliced read
        bool spliced;
        
        // Location of the alignment relative to the chromosome
        Locus l;

        // Only valid if the alignment is spliced
        Base skipped;
    };    
}

#endif