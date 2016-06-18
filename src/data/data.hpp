#ifndef DATA_HPP
#define DATA_HPP

#include "data/locus.hpp"
#include "data/types.hpp"

namespace Anaquin
{
    struct TransData_
    {
        // Eg: chr1
        ChrID cID;
        
        // Eg: ENSG00000223972.5
        GeneID gID;
        
        // Eg: ENST00000456328.2
        TransID tID;
        
        Locus l;
    };
    
    struct GeneData
    {
        // Eg: chr1
        ChrID cID;
        
        // Eg: ENSG00000223972.5
        GeneID gID;
        
        Locus l;
    };

    /*
     * Represents a mathced element that can be identified
     */
    
    struct Matched
    {
        virtual std::string name() const = 0;
    };
}

#endif