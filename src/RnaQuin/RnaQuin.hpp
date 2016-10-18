#ifndef R_RNAQUIN_HPP
#define R_RNAQUIN_HPP

#include "data/data.hpp"

namespace Anaquin
{
    static bool isRnaQuin(const ChrID &id)
    {
        return id == "chrIS";
    }
    
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
