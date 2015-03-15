#ifndef AS_STANDARD_FACTORY_HPP
#define AS_STANDARD_FACTORY_HPP

#include <map>
#include <memory>
#include "Feature.hpp"
#include "Sequence.hpp"

typedef std::map<Position, Feature> FeatureMap;

struct StandardFactory
{
	static std::string transGTF();

	// Returns name of the in-sillico chromosome
	static std::string chromoName();

    // Returns the features for the chromosome indexed with a mapping
    static FeatureMap features();
    
	// Returns a in-sillico sequence
	static std::shared_ptr<Sequence> sequence();
};

#endif