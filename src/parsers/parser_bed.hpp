#ifndef GI_PARSER_BED_HPP
#define GI_PARSER_BED_HPP

#include <vector>
#include "locus.hpp"
#include <functional>

namespace Spike
{
    struct BedFeature
    {
        ChromoID id;
        
        Locus l;
        
        FeatureName name;
        
        // Locations of the sorted blocks
        std::vector<Locus> blocks;
    };
    
    struct ParserBED
    {
        static bool parse(const std::string &file, std::function<void(const BedFeature &)>);
    };
}

#endif