#ifndef GI_SENSITIVITY_HPP
#define GI_SENSITIVITY_HPP

#include "types.hpp"

namespace Spike
{
    struct Sensitivity
    {
        SequinID id;
        Counts counts;

        // The lowest abundance that is still detectable
        Concentration abund;
    };
}

#endif