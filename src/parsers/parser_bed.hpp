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
        typedef std::string Feature;
        
        struct Data
        {
            operator const Feature &() const { return name; }
            
            ChromoID id;
            
            // Forward or reverse strand?
            Strand strand;

            Locus l;

            Feature name;
        };

        typedef std::function<void(const Data &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif