#ifndef GI_FEATURE_HPP
#define GI_FEATURE_HPP

#include "locus.hpp"

namespace Spike
{
    enum FeatureType
    {
        CDS,
        Exon,
        Junction,
        StopCodon,
        StartCodon,
        Transcript,
    };

    typedef std::string FeatureID;
    
    struct Feature
    {
        inline bool overlap(const Locus &l) const
        {
            return this->l.overlap(l);
        }

        operator Locus &()
        {
            return l;
        }

        void operator=(const Feature &f)
        {
            l  = f.l;
            id = f.id;
            type = f.type;
            iID  = f.iID;
            geneID = f.geneID;
        }
        
        FeatureID id;

        // The location of the feature relative to the chromosome
        Locus l;
        
        FeatureType type;
        
        // Empty if the information is unavailable
        GeneID geneID;

        // Empty if the information is unavailable
        TranscriptID iID;
    };
}

#endif