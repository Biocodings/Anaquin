#ifndef PARSER_CDIFFS_HPP
#define PARSER_CDIFFS_HPP

#include "data/dtest.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserCDiffs
    {
        typedef std::string TrackID;
        
        struct Data : public DiffTest
        {
            TrackID tID;
            
            // Test statistics
            double stats;
        };

        static void parse(const FileName &, std::function<void (const Data &, const ParserProgress &)>);
    };    
}

#endif
