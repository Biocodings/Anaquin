#ifndef GI_M_BLAST_HPP
#define GI_M_BLAST_HPP

#include <set>
#include "data/locus.hpp"
#include "data/types.hpp"
#include "data/sequin.hpp"

namespace Spike
{
    /*
     * Represents all alignments for a particular metaquin
     */

    struct MetaAlignment
    {
        // Name of the metaquin
        std::string id;

        // Mixture A and B
        Sequin seqA, seqB;

        // Contigs aligned to this metaquin
        std::vector<std::string> ids;
        
        // Contigs aligned to this metaquin
        std::vector<Locus> aligns;

        std::vector<std::string> temp;
        
        // Fraction of bases covered by alignments
        double coverage;
        
        // Fraction of bases not covered by alignments
        double mismatch;

        // Fraction of gap bases in alignments
        double gaps;
        
        inline bool operator<(const MetaAlignment &x)  const { return id < x.id;  }
        inline bool operator==(const MetaAlignment &x) const { return id == x.id; }
    };

    struct MBlast
    {
        struct Stats
        {
            std::set<MetaAlignment> metas;
            std::map<ContigID, MetaAlignment> aligns;
        };

        static Stats analyze(const std::string &file);
    };
}

#endif