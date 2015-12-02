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
                
                std::map<ExonID,   Counts> eContains, eOverlaps;
                std::map<IntronID, Counts> iContains, iOverlaps;

                // Number of exons detected for each gene
                std::map<GeneID, Counts> detectExons;
                
                // Number of exons undetected for each gene
                std::map<GeneID, Counts> undetectExons;

                // Number of introns detected for each gene
                std::map<GeneID, Counts> detectIntrons;

                // Number of introns undetected for each gene
                std::map<GeneID, Counts> undetectIntrons;

                Counts eMapped  = 0;
                Counts iMapped  = 0;
                Counts eUnknown = 0;
                Counts iUnknown = 0;

                // Intervals for exons in TransQuin
                Intervals<TransRef::ExonInterval> eInters;
                
                // Intervals for introns in TransQuin
                Intervals<TransRef::IntronInterval> iInters;

                // Performance for the alignments at exon and intron level
                Confusion ae, ai;
                
                // Overall performance at various levels
                Performance pb, pe, pi;

                // Individual performance for each sequin
                std::map<GeneID, Confusion> sb, se, si;

                // Sequins that have failed to be detected
                std::vector<MissingSequin> missings;

                // Alignments that have no mapping
                std::vector<UnknownAlignment> unknowns;

                // Overall sensitivity
                inline double sn(enum AlignMetrics l) const
                {
                    switch (l)
                    {
                        case Exon:   { return pe.m.sn(); }
                        case Intron: { return pi.m.sn(); }
                        case Base:   { return pb.m.sn(); }
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
                        case Exon:   { return ae.ac(); }
                        case Intron: { return ai.ac(); }
                        case Base:   { return pb.m.ac(); }
                    }
                }

                // Sensitivity at the gene level
                inline double sn(const GeneID &id) const
                {
                    /*
                     * Sensitivity for the gene level is defined as:
                     *
                     *                 detected exons
                     *     ------------------------------------
                     *        detected exons + undetected exons
                     *
                     * For example, if we're able to detect half of the exons in a gene.
                     * The sensitivty would be 50%.
                     */
                    
                    const auto succeed = detectExons.at(id);
                    const auto failed  = undetectExons.at(id);
                    
                    return static_cast<double>(succeed) / (succeed + failed);
                }
            };

            static Stats stats (const FileName &, const Options &options = Options());
            static Stats stats (const std::vector<Alignment> &, const Options &options = Options());

            static void  report(const FileName &, const Options &options = Options());
    };
}

#endif