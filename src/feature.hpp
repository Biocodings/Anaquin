#ifndef GI_FEATURE_HPP
#define GI_FEATURE_HPP

#include <map>
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
    
    typedef std::string OptionID;
    typedef std::string FeatureID;
    typedef std::string OptionValue;
    
    // Eg: "exons": "my_exon_name"
    typedef std::map<OptionID, OptionValue> Options;
    
    struct Feature
    {
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