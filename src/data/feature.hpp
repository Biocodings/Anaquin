#ifndef FEATURE_HPP
#define FEATURE_HPP

#include "data/locus.hpp"
#include "data/biology.hpp"

namespace Anaquin
{
    struct Feature
    {
        inline bool overlap(const Locus &l) const
        {
            return this->l.overlap(l);
        }

        operator Locus() const { return l; }

        void operator=(const Feature &f)
        {
            l    = f.l;
            type = f.type;
            tID  = f.tID;
            gID  = f.gID;
        }
        
        ChrID cID;

        // Forward or reverse strand?
        Strand strand;
        
        // The location of the feature relative to the chromosome
        Locus l;

        RNAFeature type;
        
        // Empty if the information is unavailable
        GeneID gID;

        // Empty if the information is unavailable
        TransID tID;
    };
}

#endif