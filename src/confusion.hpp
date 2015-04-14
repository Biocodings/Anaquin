#ifndef GI_CONFUSION_HPP
#define GI_CONFUSION_HPP

#include <math.h>
#include "types.hpp"

namespace Spike
{
    struct Confusion
    {
        Percentage fp = 0;
        Percentage tp = 0;
        Percentage fn = 0;
        Percentage tn = 0;
        
        inline Percentage sp() const
        {
            return (tp + fn) ? tp / (tp + fn) : NAN;
        }
        
        inline Percentage sn() const
        {
            return (fp + tn) ? tn / (fp + tn) : NAN;
        }
    };    
}

#endif