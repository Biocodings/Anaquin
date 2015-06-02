#ifndef GI_PARSER_BLAT_HPP
#define GI_PARSER_BLAT_HPP

#include "types.hpp"
#include <functional>
#include "parsers/parser.hpp"

namespace Spike
{
    struct ParserBlat
    {
        struct BlatLine
        {
            // Target sequence name
            std::string tName;

            // Query sequence name
            std::string qName;
            
            // Alignment start position in query
            BasePair tStart;

            // Alignment end position in query
            BasePair tEnd;
        };

        typedef std::function<void(const BlatLine &, const ParserProgress &)> Callback;
        static void parse(const std::string &, Callback, ParserMode mode = File);
    };
}

#endif