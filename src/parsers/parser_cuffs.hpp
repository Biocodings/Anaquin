#ifndef GI_PARSER_CUFFS_HPP
#define GI_PARSER_CUFFS_HPP

#include <map>
#include "data/types.hpp"

namespace Anaquin
{
    enum TrackingStatus
    {
        OK,
        HIData,
        NoTest
    };

    typedef std::string TrackID;

    static const std::map<TrackID, TrackingStatus> tok2Status =
    {
        { "OK", OK },
        { "HIDATA", HIData },
        { "NOTEST", NoTest }
    };
}

#endif
