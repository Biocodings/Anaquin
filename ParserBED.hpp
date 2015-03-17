#ifndef AS_PARSER_BED_HPP
#define AS_PARSER_BED_HPP

#include <vector>
#include <functional>
#include "Feature.hpp"

struct BedFeature : public Feature
{
    std::string name;

    // Locations of the sorted blocks
    std::vector<Locus> blocks;
};

struct ParserBED
{
    static bool parse(const std::string &file, std::function<void(const BedFeature &)>);
};

#endif