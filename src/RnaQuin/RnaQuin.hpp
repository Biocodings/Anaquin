#ifndef R_RNAQUIN_HPP
#define R_RNAQUIN_HPP

#include "data/data.hpp"

namespace Anaquin
{
    const ChrID ChrIS = "chrIS";

    inline bool isChrIS(const ChrID &x) { return x == ChrIS; }

    // Eg: R1_1_1 to R1_1
    inline SequinID isoform2Gene(const SequinID &x)
    {
        return x.substr(0, x.size()-2);
    }
}

#endif
