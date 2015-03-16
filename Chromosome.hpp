#ifndef AS_CHROMOSOME_HPP
#define AS_CHROMOSOME_HPP

#include <map>
#include <string>
#include <vector>
#include "Feature.hpp"

typedef std::map<std::string, Sequence> SequinMap;

struct Chromosome
{
    std::string id;
    std::vector<Feature> fs;
    
    // The list of sequins added to the experiment
    SequinMap sequins;
};

#endif