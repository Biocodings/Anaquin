#ifndef PARSER_BLAST_HPP
#define PARSER_BLAST_HPP

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
            
            // Alignment start position in target
            Base tStart;

            // Alignment end position in target
            Base tEnd;
            
            // Target sequence size
            Base tSize;
            
            // Alignment start position in query
            Base qStart;
            
            // Alignment end position in query
            Base qEnd;
            
            // Query sequence size
            Base qSize;

            // Number of inserts in query
            Counts qGapCount;
            
            // Number of bases inserted into query
            Base qGap;
            
            // Number of inserts in target
            Counts tGapCount;
            
            // Number of bases inserted into target
            Base tGap;
            
            // Number of matching bases
            Base match;

            // Number of mismatching bases
            Base mismatch;
        };

        typedef std::function<void(const BlastLine &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif