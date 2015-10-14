#ifndef PARSER_SAM_HPP
#define PARSER_SAM_HPP

#include "data/alignment.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserSAM
    {
        struct AlignmentInfo
        {
            ParserProgress p;

            // Size of the chromosome of the alignment
            Base size;
        };
        
        typedef std::function<void (const Alignment &, const AlignmentInfo &)> Callback;
        
        static void parse(const FileName &file, std::function<void (const Alignment &, const AlignmentInfo &)>);
    };
}

#endif