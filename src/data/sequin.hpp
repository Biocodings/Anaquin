#ifndef GI_SEQUIN_HPP
#define GI_SEQUIN_HPP

#include "data/types.hpp"
#include "data/locus.hpp"

namespace Spike
{
    struct Sequin
    {
        operator Locus()     const { return l;  }
        operator IsoformID() const { return id; }
        
        inline Concentration abund() const { return raw; }
        
        SequinID id;
        Locus l;
        
        // Amount of abundance
        Concentration raw;
    };
    
    struct Sequins
    {
        inline Concentration abund() const
        {
            return r.abund() + v.abund();
        }
        
        // Each mixture represents a transcript for a gene
        GeneID geneID;
        
        // Reference and variant mixtures
        Sequin r, v;
    };
}

#endif