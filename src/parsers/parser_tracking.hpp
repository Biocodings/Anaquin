#ifndef PARSER_TRACKING_HPP
#define PARSER_TRACKING_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cuffs.hpp"

namespace Anaquin
{
    struct ParserTracking
    {
        struct Data : public Expression
        {
            TrackID trackID;
            
            FPKM lFPKM;
            FPKM uFPKM;
            
            TrackingStatus status;
        };

        static void parse(const FileName &, std::function<void (const Data &, const ParserProgress &)>);
    };    
}

#endif