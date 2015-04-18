#ifndef GI_CONFUSION_HPP
#define GI_CONFUSION_HPP

#include <math.h>
#include "types.hpp"

namespace Spike
{
    struct Confusion
    {
        Counts fp = 0;
        Counts tp = 0;
        Counts fn = 0;

        inline Percentage sp() const
        {
            return (tp + fp) ? tp / (tp + fp) : NAN;
        }

        inline Percentage sn() const
        {
            return (tp + fn) ? tp / (tp + fn) : NAN;
        }
    };    
}

#endif