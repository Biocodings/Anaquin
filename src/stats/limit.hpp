#ifndef LIMIT_HPP
#define LIMIT_HPP

#include <math.h>
#include "data/data.hpp"

namespace Anaquin
{
    struct Limit
    {
        SequinID id;

        // Expected concentration for the sequin
        Concent abund = NAN;
    };
}

#endif