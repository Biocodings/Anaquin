#ifndef GI_PARSER_CUFFS_HPP
#define GI_PARSER_CUFFS_HPP

#include <map>
#include "types.hpp"

namespace Spike
{
    enum TrackingStatus
    {
        OK,
        HIData
    };

    typedef std::string TrackID;
    
    static const std::map<TrackID, TrackingStatus> tok2Status =
    {
        { "OK", OK },
        { "HIDATA", HIData }
    };
}

#endif
