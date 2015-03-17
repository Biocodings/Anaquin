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

typedef std::map<std::string, std::string> Options;

struct Feature
{
	std::string id;

	Locus end;
	Locus start;
	Locus length;

    FeatureType type;

    // Optional field such as "gene_id" and "transcript_id"
    Options options;
};

#endif