#ifndef PARSER_BED_HPP
#define PARSER_BED_HPP

#include <functional>
#include "data/types.hpp"
#include "data/locus.hpp"
#include "data/reader.hpp"
#include "data/biology.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserBed
    {
        struct Data
        {
            operator const FeatureName &() const { return name; }
            
            ChromoID id;
            
            // Forward or reverse strand?
            Strand strand;

            Locus l;
            
            FeatureName name;
            
            // Locations of the sorted blocks
            //std::vector<Locus> blocks;
        };

        typedef std::function<void(const Data &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif