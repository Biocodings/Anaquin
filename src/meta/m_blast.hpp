#ifndef GI_M_BLAST_HPP
#define GI_M_BLAST_HPP

#include <set>
#include "data/types.hpp"

namespace Spike
{
    struct MetaAlignment
    {
        // Name of the metaquin
        std::string id;
    };
    
    struct MBlast
    {
        struct Stats
        {
            std::set<MetaAlignment> aligns;
        };

        static Stats analyze(const std::string &file);
    };
}

#endif