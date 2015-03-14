#ifndef AS_SILLICO_FACTORY_HPP
#define AS_SILLICO_FACTORY_HPP

#include <map>
#include <memory>
#include "Feature.hpp"
#include "Sequence.hpp"

typedef std::map<Position, Feature> FeatureMap;

/*
 * This factory class provides support for the in-sillico chromosome developed by
 * Gavian Institute of Medical Research.
 */

struct SillicoFactory
{
	static std::string transGTF();

    // Returns the features for the chromosome indexed with a mapping
    static FeatureMap features();
    
	// Returns a in-sillico sequence
	static std::shared_ptr<Sequence> sequence();
};

#endif