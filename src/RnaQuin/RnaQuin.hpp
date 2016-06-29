#ifndef R_RNAQUIN_HPP
#define R_RNAQUIN_HPP

#include "data/data.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RnaQuin
    {
        // Eg: R1_1_1 to R1_1
        static SequinID t2g(const SequinID &id)
        {
            return id.substr(0, id.size()-2);
        }
    };
}

#endif