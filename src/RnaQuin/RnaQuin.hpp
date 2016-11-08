#ifndef R_RNAQUIN_HPP
#define R_RNAQUIN_HPP

#include "data/data.hpp"

namespace Anaquin
{
    inline bool isRnaQuin(const ChrID &id)
    {
        return id == ChrIS;
    }

    // Eg: R1_1_1 to R1_1
    inline SequinID Isoform2Gene(const SequinID &id)
    {
        return id.substr(0, id.size()-2);
    }
}

#endif
