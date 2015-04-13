#ifndef GI_PARSER_TRACKING_HPP
#define GI_PARSER_TRACKING_HPP

#include "types.hpp"
#include <functional>

namespace Spike
{
    enum CTrackingStatus
    {
        OK,
        HIData
    };

    typedef std::string TrackID;
    
    struct Tracking
    {
        TrackID trackID;
        TrackID geneID;
        
        FPKM fpkm;
        FPKM lFPKM;
        FPKM uFPKM;
        
        CTrackingStatus status;
    };
    
    struct ParserCTracking
    {
        static bool parse(const std::string &file, std::function<void (const Tracking &)>);
    };    
}

#endif
