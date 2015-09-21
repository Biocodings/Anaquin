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
        Base match;

        // Number of gaps in the sequin
        Base gap;

        // Number of mis-matching bases
        Base mismatch;
    };

    /*
     * Represents alignments for a particular sequin
     */

    struct MetaAlignment
    {
        inline const SequinID &id() const { return seq->id; }
        
        const SequinData *seq;

        // Contigs aligned to this metaquin
        std::vector<AlignedContig> contigs;

        /*
         * The following metrics are only valid if an alignment is available
         */

        // Proportion of non-overlapping bases covered by alignments
        double covered = 0.0;
        
        // Proportion of bases not covered by alignments
        double oMismatch = 0.0;

        // Proportion of gap bases in alignments
        double oGaps = 0.0;

        // Average coverage depth across assembled sequence
        double depthAlign = 0.0; // TODO: ????

        // Average coverage depth across entire sequin
        double depthSequin = 0.0; // TODO: ????

        inline bool operator<(const MetaAlignment &x)  const { return seq->id < x.seq->id;  }
        inline bool operator==(const MetaAlignment &x) const { return seq->id == x.seq->id; }
    };

    typedef std::map<SequinID, std::shared_ptr<MetaAlignment>> SequinAlign;

    struct MBlat
    {
        struct Stats : public MappingStats
        {
            // Proportion of overlapping bases assembled for all sequins
            inline double overMatch() const { return static_cast<double>(oMatch) / total; }

            // Proportion of overlapping gaps for all sequins
            inline double overGaps() const { return static_cast<double>(oGaps) / total; }

            // Proportion of overlapping mismatches for all sequins
            inline double overMismatch() const { return static_cast<double>(oMismatch) / total; }

            // Total number of sequins that are assembled (i.e: at least a single contig)
            inline Counts countAssembled() const
            {
                return std::count_if(metas.begin(),metas.end(),
                            [&](const std::pair<SequinID, std::shared_ptr<MetaAlignment>> &p)
                            {
                                return p.second->contigs.empty() ? 0 : 1;
                            });
            }

            // Sum of bases for all sequins
            Base total = 0;

            // Sum of overlapping gaps for all sequins
            Base oGaps = 0;

            // Sum of overlapping mismatch for all sequins
            Base oMismatch = 0;
            
            // Sum of overlapping matching bases for all sequins
            Base oMatch = 0;
            
            // For each sequin (could be unmapped)
            std::map<SequinID, std::shared_ptr<MetaAlignment>> metas;

            // For each contig listed in the alignment file
            std::map<ContigID, std::shared_ptr<MetaAlignment>> aligns;
        };

        typedef AnalyzerOptions Options;
        
        /*
         * Conduct statistical analysis for an alignment geneated by BLAST against
         * the sequins. The input file is assumed to be PSL format.
         */

        static Stats analyze(const FileName &, const Options &o = Options());

        /*
         * Generate summary statistics for the alignment
         */

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif