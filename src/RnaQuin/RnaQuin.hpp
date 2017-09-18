#ifndef RNAQUIN_HPP
#define RNAQUIN_HPP

#include "data/data.hpp"

namespace Anaquin
{
    inline ChrID ChrIS() { return "chrIS"; }

    inline bool isChrIS(const ChrID &x) { return x == "chrIS" || x == "IS"; }

    // Eg: R1_1_1 to R1_1
    inline SequinID isoform2Gene(const SequinID &x)
    {
        return x.substr(0, x.size()-2);
    }
}

#endif
