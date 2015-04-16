#ifndef GI_SENSITIVITY_HPP
#define GI_SENSITIVITY_HPP

#include "types.hpp"
#include "expression.hpp"

namespace Spike
{
    struct Sensitivity
    {
        Sensitivity() {}
        Sensitivity(SequinID id, Counts counts, Concentration exp)
                : id(id), counts(counts), exp(exp) {}

        SequinID id;
        Counts counts;
        Concentration exp;
    };
}

#endif