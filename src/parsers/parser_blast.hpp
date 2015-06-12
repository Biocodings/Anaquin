#ifndef GI_PARSER_BLAST_HPP
#define GI_PARSER_BLAST_HPP

#include <functional>
#include "data/types.hpp"
#include "parsers/parser.hpp"

namespace Spike
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
            BasePair tStart;

            // Alignment end position in query
            BasePair tEnd;

            // Number of gap bases in query
            BasePair qGaps;
            
            // Number of gap bases in target
            BasePair tGaps;
            
            // Number of matching bases
            BasePair matches;

            // Number of mismatching bases
            BasePair mismatch;
        };

        typedef std::function<void(const BlastLine &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif