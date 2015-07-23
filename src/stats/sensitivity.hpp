#ifndef GI_SENSITIVITY_HPP
#define GI_SENSITIVITY_HPP

#include "data/types.hpp"

namespace Anaquin
{
    struct Sensitivity
    {
        SequinID id;

        // Measured abundance for the limited sequin
        Counts counts;

        // Expected abundance for the limited sequin
        Concentration abund;
    };
}

#endif