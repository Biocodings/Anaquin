#ifndef AS_STANDARD_FACTORY_HPP
#define AS_STANDARD_FACTORY_HPP

#include <map>
#include "Sequence.hpp"
#include "Chromosome.hpp"

typedef std::map<std::string, Sequence> SequinMap;

struct StandardFactory
{
    // Returns the in-sillico chromosome constructed from sequins
    static Chromosome reference();

    // Returns list of sequins mixed with a sample
    static SequinMap sequins();
};

#endif