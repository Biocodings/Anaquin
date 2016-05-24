#ifndef PARSER_SAM_HPP
#define PARSER_SAM_HPP

#include "data/alignment.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserSAM
    {
        struct Info
        {
            ParserProgress p;

            // Whether this is a spliced-alignment
            bool spliced;
            
            // Size of the chromosome of the alignment
            Base length;
            
            // Internal data representation
            void *data;
            
            // Internal data representation
            void *header;
        };
        
        typedef std::function<void (const Alignment &, const Info &)> Functor;

        static void parse(const FileName &file, Functor);
    };
}

#endif