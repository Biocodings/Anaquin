#ifndef GI_M_BLAST_HPP
#define GI_M_BLAST_HPP

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
        // For convenience, the outputs are done in either direction
        struct Stats
        {
            // For each contig in the mixture (could be unmapped)
            std::map<SequinID, MetaAlignment> metas;

            // For each contig listed in the alignment file
            std::map<ContigID, MetaAlignment> aligns;
        };

        /*
         * Conduct statistical analysis for an alignment geneated by BLAST relative
         * the sequins. The input file is assumed to be PSL format.
         */

        static Stats analyze(const std::string &);
    };
}

#endif