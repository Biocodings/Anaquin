#ifndef GI_SENSITIVITY_HPP
#define GI_SENSITIVITY_HPP

#include "types.hpp"
#include "standard.hpp"
#include "expression.hpp"

namespace Spike
{
    struct Sensitivity
    {
        Sensitivity() {}
        Sensitivity(const Standard &s, const Expression::ExpressionResults<IsoformID> &r)
        {
            id = r.limit_key;
            counts = r.limit_count;
            exp = r.limit_count ? s.seqs_iA.at(r.limit_key).exp +
                                  s.seqs_iA.at(r.limit_key).exp: NAN;
        }

        SequinID id;
        Counts counts;
        Concentration exp;
    };
}

#endif