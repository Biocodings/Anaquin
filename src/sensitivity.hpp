#ifndef GI_SENSITIVITY_HPP
#define GI_SENSITIVITY_HPP

#include "types.hpp"

namespace Spike
{
    struct Sensitivity
    {
        SequinID id;
        
        Counts counts;
        
        Concentration abundance;
    };
}

#endif