#ifndef PARSER_TRACKING_HPP
#define PARSER_TRACKING_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cuffs.hpp"

namespace Anaquin
{
    struct Tracking : public Expression
    {
        TrackID trackID;

        FPKM lFPKM;
        FPKM uFPKM;
        
        TrackingStatus status;
    };
    
    struct ParserTracking
    {
        static void parse(const std::string &, std::function<void (const Tracking &, const ParserProgress &)>);
    };    
}

#endif