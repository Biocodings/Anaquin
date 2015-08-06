#ifndef GI_PARSER_BLAST_HPP
#define GI_PARSER_BLAST_HPP

#include <functional>
#include "data/types.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserBlast
    {
        struct BlastLine
        {
            // Target sequence name
            std::string tName;

            // Query sequence name
            std::string qName;
            
            // Alignment start position in query
            Base tStart;

            // Alignment end position in query
            Base tEnd;

            // Number of gap bases in query
            Base qGaps;
            
            // Number of gap bases in target
            Base tGaps;
            
            // Number of matching bases
            Base matches;

            // Number of mismatching bases
            Base mismatch;
        };

        typedef std::function<void(const BlastLine &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif