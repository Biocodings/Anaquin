#ifndef PARSER_CDIFFS_HPP
#define PARSER_CDIFFS_HPP

#include "data/dtest.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser_cuffs.hpp"

namespace Anaquin
{
    struct TrackingDiffs : public DiffTest
    {
        TrackID testID;

        // Test statistics
        double stats;
    };

    struct ParserCDiffs
    {
        static void parse(const FileName &, std::function<void (const TrackingDiffs &, const ParserProgress &)>);
    };    
}

#endif
