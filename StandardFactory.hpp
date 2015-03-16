#ifndef AS_STANDARD_FACTORY_HPP
#define AS_STANDARD_FACTORY_HPP

#include <map>
#include "Feature.hpp"
#include "Sequence.hpp"
#include "Chromosome.hpp"

struct StandardFactory
{
    static Chromosome reference();

	// Returns a in-sillico sequence
	static std::shared_ptr<Sequence> sequence();
};

#endif