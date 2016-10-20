#ifndef M_BLAT_HPP
#define M_BLAT_HPP

#include "stats/analyzer.hpp"
#include "data/reference.hpp"

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
        
        // Number of mismatching bases
        Base mismatch;
        
        // Number of gaps in the sequin
        Base rGap;
        
        // Alignment start position in sequin
        Base rStart;
        
        // Alignment end position in sequin
        Base rEnd;
        
        // Target sequin size
        Base rSize;
        
        // Alignment start position in query
        Base qStart;
        
        // Alignment end position in query
        Base qEnd;
        
        // Number of gaps in the query
        Base qGap;
        
        // Query sequence size
        Base qSize;
        
        // Number of inserts in query
        Counts qGapCount;
        
        // Number of inserts in target
        Counts rGapCount;
    };
    
    /*
     * Represents contig alignments for a sequin
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
        Proportion covered = 0.0;
        
        Proportion oMatch = 0.0;
        
        // Proportion of bases not covered by alignments
        Proportion oMismatch = 0.0;
        
        // Proportion of bases that are gaps in the sequin (reference)
        Proportion oRGaps = 0.0;
        
        // Proportion of bases that are gaps in all contigs aligned to this sequin
        Proportion oQGaps = 0.0;
        
        // Average coverage depth across assembled sequence
        Proportion depthAlign = 0.0; // TODO: ????
        
        // Average coverage depth across entire sequin
        Proportion depthSequin = 0.0; // TODO: ????
        
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
            
            // Proportion of overlapping gaps for all sequins (reference)
            inline double overRGaps() const { return static_cast<double>(oRGaps) / total; }
            
            // Proportion of overlapping gaps for all sequins (reference)
            inline double overQGaps() const { return static_cast<double>(oQGaps) / oQSums; }
            
            // Proportion of overlapping mismatches for all sequins
            inline double overMismatch() const { return static_cast<double>(oMismatch) / total; }
            
            // Total number of sequins that are assembled (i.e: at least a single contig)
            inline Counts countAssembled() const
            {
                return std::count_if(metas.begin(), metas.end(), [&](const std::pair<SequinID, std::shared_ptr<MetaAlignment>> &p)
                                     {
                                         return p.second->contigs.empty() ? 0 : 1;
                                     });
            }
            
            // Sum of bases for all sequins
            Base total = 0;
            
            // Sum of overlapping gaps for all sequins (ie: references)
            Base oRGaps = 0;
            
            // Sum of overlapping gaps for all contigs aligned
            Base oQGaps = 0;
            
            // Sum of overlapping sizes for all contigs aligned
            Base oQSums = 0;
            
            // Sum of overlapping mismatch for all sequins
            Base oMismatch = 0;
            
            // Sum of overlapping matching bases for all sequins
            Base oMatch = 0;
            
            // For each sequin (could be unmapped)
            std::map<SequinID, std::shared_ptr<MetaAlignment>> metas;
            
            // For each contig listed in the alignment file
            std::map<ContigID, std::shared_ptr<MetaAlignment>> aligns;
            
            /*
             * Mappings for targets (including sequins)
             */
            
            std::map<SequinID, Base> t2l;
            
            /*
             * Mapping for contigs
             */
            
            // Mapping from sequin to contig
            std::map<ContigID, SequinID> c2s;
            
            // Mapping from contig to their aligned (target) length
            std::map<ContigID, Base> c2a;
            
            // Mapping from contig to their length (size of the entire contig)
            std::map<ContigID, Base> c2l;
        };
        
        typedef AnalyzerOptions Options;
        
        /*
         * Conduct statistical analysis for an alignment geneated by BLAST against
         * the sequins. The input file is assumed to be PSL format.
         */
        
        static Stats analyze(const FileName &, const Options &o = Options());
    };
}

#endif
