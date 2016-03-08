#ifndef LIMIT_HPP
#define LIMIT_HPP

#include <math.h>
#include "data/types.hpp"

namespace Anaquin
{
    struct Limit
    {
        SequinID id;

        // Measured abundance for the limited sequin
        Counts counts;

        // Expected abundance for the limited sequin
        Concent abund = NAN;
    };
}

#endif