#ifndef GI_PARSER_BED_HPP
#define GI_PARSER_BED_HPP

#include <functional>
#include "data/types.hpp"
#include "data/locus.hpp"
#include "parsers/parser.hpp"

namespace Spike
{
    struct BedFeature
    {
        operator const FeatureName &() const { return name; }

        ChromoID id;
        
        Locus l;

        FeatureName name;
        
        // Locations of the sorted blocks
        std::vector<Locus> blocks;
    };

    struct ParserBED
    {
        typedef std::function<void(const BedFeature &, const ParserProgress &)> Callback;
        static void parse(const std::string &, Callback, DataMode mode = File);
    };
}

#endif