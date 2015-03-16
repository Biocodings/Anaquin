#ifndef AS_STANDARD_FACTORY_HPP
#define AS_STANDARD_FACTORY_HPP

#include <map>
#include "Sequence.hpp"
#include "Chromosome.hpp"

typedef std::map<std::string, Sequence> SequinMap;

struct StandardFactory
{
    static Chromosome reference();

    static SequinMap sequins();

	// Returns a in-sillico sequence
	static std::shared_ptr<Sequence> sequence();
};

#endif