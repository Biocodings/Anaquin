#ifndef GI_PARSER_BED_HPP
#define GI_PARSER_BED_HPP

#include "types.hpp"
#include "locus.hpp"
#include <functional>
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
        static void parse(const std::string &, Callback, ParserMode mode = File);
    };
}

#endif