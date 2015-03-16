#ifndef AS_CHROMOSOME_HPP
#define AS_CHROMOSOME_HPP

#include <vector>
#include "Feature.hpp"

/*
 * This class abstracts information required by a chromosome. However, the sequence is not included.
 */

struct Chromosome
{
    std::string id;
    std::vector<Feature> fs;
};

#endif