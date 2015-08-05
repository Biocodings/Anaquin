#ifndef GI_M_PSL_HPP
#define GI_M_PSL_HPP

#include "data/sequin.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    /*
     * Represents a contig that has been aligned.
     */

    struct AlignedContig
    {
        operator const Locus &() const { return l; }

        ContigID id;

        // The positon where the alignment occurs
        Locus l;

        // Number of matching bases
        BasePair match;

        // Number of gaps in the sequin
        BasePair gap;

        // Number of mis-matching bases
        BasePair mismatch;
    };

    /*
     * Represents alignments for a particular sequin
     */

    struct MetaAlignment
    {
        // Name of the metaquin
        SequinID id;

        // Mixture A and B
        Sequin seqA, seqB;

        // Contigs aligned to this metaquin
        std::vector<AlignedContig> contigs;

        /*
         * The following metrics are only valid if an alignment is available
         */

        // Fraction of bases covered by alignments
        double covered = 0.0;
        
        // Fraction of bases not covered by alignments
        double mismatch = 0.0;

        // Fraction of gap bases in alignments
        double gaps = 0.0;

        // Average coverage depth across assembled sequence
        double depthAlign = 0.0; // TODO: ????
        
        // Average coverage depth across entire sequin
        double depthSequin = 0.0; // TODO: ????

        inline bool operator<(const MetaAlignment &x)  const { return id < x.id;  }
        inline bool operator==(const MetaAlignment &x) const { return id == x.id; }
    };

    typedef std::map<SequinID, MetaAlignment> SequinAlign;

    struct MBlast
    {
        // For convenience, the outputs are done in either direction
        struct Stats
        {
            // For each sequin (could be unmapped)
            std::map<SequinID, MetaAlignment> metas;

            // For each contig listed in the alignment file
            std::map<ContigID, MetaAlignment> aligns;
        };

        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        /*
         * Conduct statistical analysis for an alignment geneated by BLAST relative
         * the sequins. The input file is assumed to be PSL format.
         */

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif