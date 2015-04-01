#ifndef GI_FEATURE_HPP
#define GI_FEATURE_HPP

#include <map>
#include "locus.hpp"

enum FeatureType
{
	CDS,
    Exon,
    Intron,
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

    // Empty if the information is not available
    GeneID geneID;
    
    // Empty if the information is not available
    IsoformID iID;
    
    // Optional field such as "gene_id" and "transcript_id"
    //Options options;
};

#endif
