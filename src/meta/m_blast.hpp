#ifndef GI_M_BLAST_HPP
#define GI_M_BLAST_HPP

#include <set>
#include "data/locus.hpp"
#include "data/types.hpp"

namespace Spike
{
    struct MetaAlignment
    {
        // Name of the metaquin
        std::string id;

        Concentration mixA, mixB;

        // Contigs aligned to this metaquin
        std::set<Locus> aligns;

        inline bool operator==(const MetaAlignment &x) const { return id == x.id; }
        inline bool operator<(const MetaAlignment &x)  const { return id < x.id;  }
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