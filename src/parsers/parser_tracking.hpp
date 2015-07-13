#ifndef GI_PARSER_TRACKING_HPP
#define GI_PARSER_TRACKING_HPP

#include "stats/analyzer.hpp"
#include "data/types.hpp"
#include "parsers/parser_cuffs.hpp"

namespace Anaquin
{
    typedef std::string TrackID;
    
    struct Tracking
    {
        TrackID trackID;
        TrackID geneID;

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
