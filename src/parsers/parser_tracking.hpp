#ifndef GI_PARSER_TRACKING_HPP
#define GI_PARSER_TRACKING_HPP

#include "types.hpp"
#include <functional>

namespace Spike
{
    enum TrackingStatus
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
        
        TrackingStatus status;
    };
    
    struct ParserTracking
    {
        static bool parse(const std::string &file, std::function<void (const Tracking &)>);
    };    
}

#endif
