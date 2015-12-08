#ifndef T_ALIGN_HPP
#define T_ALIGN_HPP

#include "stats/analyzer.hpp"
#include "data/alignment.hpp"

namespace Anaquin
{
    class TAlign : public Analyzer
    {
        public:
            typedef AnalyzerOptions Options;

            typedef std::map<BinID, Counts> BinCounts;

            struct MergedConfusion
            {
                /*
                 * Alignment statistics
                 */
                
                Counts aTP = 0;
                Counts aFP = 0;
            
                /*
                 * Level statistics (eg: exons and introns)
                 */
                
                Counts lTP = 0;
                Counts lNR = 0;

                inline Counts aNQ() const { return aTP + aFP; }
                inline Counts lFN() const { return lNR - lTP; }

                inline double sn() const
                {
                    return (lTP + lFN()) ? static_cast<double>(lTP) / (lTP + lFN()) : NAN;
                }
                
                inline double ac() const
                {
                    return (aTP + aFP) ? static_cast<double>(aTP) / (aTP + aFP) : NAN;
                }
            };

            struct Stats : public AlignmentStats
            {
                enum AlignMetrics
                {
                    Exon,
                    Intron,
                    Base
                };
                
                std::map<ExonID,   GeneID> exonToGene;
                std::map<IntronID, GeneID> intronToGene;
                
                BinCounts eContains, eOverlaps;
                BinCounts iContains, iOverlaps;

                // Intervals for exons in TransQuin
                Intervals<TransRef::ExonInterval> eInters;
                
                // Intervals for introns in TransQuin
                Intervals<TransRef::IntronInterval> iInters;

                /*
                 * Overall statistics
                 */
                
                MergedConfusion overE, overI;

                Hist  histE,  histI;
                Limit limitE, limitI;
                
                /*
                 * Individual statistics for each gene (due to alternative splicing)
                 */
                
                std::map<GeneID, MergedConfusion> geneE, geneI;
                
                // Overall performance at various levels
                Performance pb;

                // Individual performance for each sequin
                std::map<GeneID, Confusion> sb;

                // Sequins that have failed to be detected
                std::vector<MissingSequin> missings;

                // Alignments that have no mapping
                std::vector<UnknownAlignment> unknowns;

                // Number of exons in the query
                inline Counts qExons() const   { return overE.aNQ();  }
                
                // Number of introns in the query
                inline Counts qIntrons() const { return overI.aNQ(); }

                // Number of bases in the query
                inline Counts qBases() const { return pb.m.nq(); }

                // Overall sensitivity
                inline double sn(enum AlignMetrics l) const
                {
                    switch (l)
                    {
                        case Exon:   { return overE.sn(); }
                        case Intron: { return overI.sn(); }
                        case Base:   { return pb.m.sn();  }
                    }
                }
                
                // Overall accuracy
                inline double ac(enum AlignMetrics l) const
                {
                    /*
                     * It's important to note pe.m.ac() and pi.m.ac() are both undefined. The two
                     * confusion matrices for used for measuring sensivitiy and measued in numbers
                     * of exons/introns.
                     */
                    
                    switch (l)
                    {
                        case Exon:   { return overE.ac(); }
                        case Intron: { return overI.ac(); }
                        case Base:   { return pb.m.ac();  }
                    }
                }

                // Sensitivity at the gene level
                inline double sn(const GeneID &id) const
                {
                    return geneE.at(id).sn();
                }
            };

            static Stats analyze(const FileName &, const Options &options = Options());
            static Stats analyze(const std::vector<Alignment> &, const Options &options = Options());
            static void  report (const FileName &, const Options &options = Options());
    };
}

#endif