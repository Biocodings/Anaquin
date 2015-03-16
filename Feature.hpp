#ifndef AS_FEATURE_HPP
#define AS_FEATURE_HPP

#include "Types.hpp"

enum FeatureType
{
	CDS,
    Exon,
    Intron,
    Junction,
	StartCodon
};

struct Feature
{
	std::string id;

	Locus end;
	Locus start;
	BasePair length;

    FeatureType type;
};

#endif