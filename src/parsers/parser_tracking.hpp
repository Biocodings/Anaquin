#ifndef GI_PARSER_TRACKING_HPP
#define GI_PARSER_TRACKING_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cuffs.hpp"

namespace Anaquin
{
    struct Tracking
    {
        TrackID chromID;
        TrackID trackID;
        TrackID geneID;

        Locus l;
        
        FPKM fpkm;
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