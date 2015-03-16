#ifndef AS_FEATURE_HPP
#define AS_FEATURE_HPP

#include "Types.hpp"

enum FeatureType
{
    Exon,
    Intron,
    Junction
};

/*
 * This class represents a biological feature, such as CDS, exon and intron.
 */

struct Feature
{
	std::string id;

    Position start;
    Position end;

    FeatureType type;
};

#endif