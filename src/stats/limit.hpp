#ifndef LIMIT_HPP
#define LIMIT_HPP

#include <math.h>
#include "data/types.hpp"

namespace Anaquin
{
    /*
     * This class represents limit of quantification
     */

    struct Limit
    {
        SequinID id;

        // Measured abundance for the sequin
        Counts counts;

        // Expected abundance for the sequin
        Concent abund = NAN;
    };
}

#endif