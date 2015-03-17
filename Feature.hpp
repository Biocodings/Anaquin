#ifndef AS_FEATURE_HPP
#define AS_FEATURE_HPP

#include <map>
#include "Types.hpp"

enum FeatureType
{
	CDS,
    Exon,
    Intron,
    Junction,
	StartCodon
};

typedef std::string FeatureID;
typedef std::string OptionID;
typedef std::string OptionValue;

// Eg: "exons": "my_exon_name"
typedef std::map<OptionID, OptionValue> Options;

struct Feature
{
	FeatureID id;

	Locus end;
	Locus start;
	Locus length;

    FeatureType type;

    // Optional field such as "gene_id" and "transcript_id"
    Options options;
};

#endif