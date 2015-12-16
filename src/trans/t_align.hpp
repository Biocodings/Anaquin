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
            typedef std::map<GeneID, Base> FPStats;

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
                
                inline double precise() const
                {
                    return (aTP + aFP) ? static_cast<double>(aTP) / (aTP + aFP) : NAN;
                }
            };

            struct Stats
            {
                enum AlignMetrics
                {
                    AlignExon,
                    AlignIntron,
                    AlignBase
                };
                
                enum MissingMetrics
                {
                    MissingExon,
                    MissingIntron,
                    MissingGene
                };
                
                struct Synthetic : public AlignmentStats
                {
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
                    
                    // Overall performance at the base level
                    Performance overB;
                    
                    Hist  histE,  histI;
                    Limit limitE, limitI;
                    
                    /*
                     * Individual statistics for each gene (due to alternative splicing)
                     */
                    
                    std::map<GeneID, MergedConfusion> geneE, geneI;
                    
                    // Individual performance at the base level
                    std::map<GeneID, Confusion> geneB;
                    
                    // Alignments that have no mapping
                    std::vector<UnknownAlignment> unknowns;
                    
                    /*
                     * Missing statistics
                     */
                    
                    std::set<Missing> missE, missI, missG;

                    // False positives (left and right)
                    FPStats lFPS, rFPS;
                };

                struct Experiment
                {
                    // Intervals for exons in TransQuin
                    Intervals<TransRef::ExonInterval> eInters;
                    
                    // Intervals for introns in TransQuin
                    Intervals<TransRef::IntronInterval> iInters;
                    
                    // Alignments that have no mapping
                    std::vector<UnknownAlignment> unknowns;
                    
                    BinCounts eContains, eOverlaps;
                    BinCounts iContains, iOverlaps;

                    /*
                     * Missing statistics
                     */
                    
                    std::set<Missing> missE, missI, missG;

                    // False positives (left and right)
                    FPStats lFPS, rFPS;
                };

                // Statistics for the synthetic chromosome
                std::shared_ptr<Synthetic> chrT;

                // Statistics for the gencode
                std::shared_ptr<Experiment> expT;
                
                /*
                 * Accessor functions
                 */

                // Number of exons in the query
                inline Counts qExons() const { return chrT->overE.aNQ(); }
                
                // Number of introns in the query
                inline Counts qIntrons() const { return chrT->overI.aNQ(); }

                // Number of bases in the query
                inline Counts qBases() const { return chrT->overB.m.nq(); }

                inline CountPercent missing(enum MissingMetrics m) const
                {
                    switch (m)
                    {
                        case MissingGene:   { return CountPercent(chrT->missG.size(), chrT->histE.size());     }
                        case MissingExon:   { return CountPercent(chrT->missE.size(), chrT->eContains.size()); }
                        case MissingIntron: { return CountPercent(chrT->missI.size(), chrT->iContains.size()); }
                    }
                }

                // Overall sensitivity
                inline double sn(enum AlignMetrics m) const
                {
                    switch (m)
                    {
                        case AlignExon:   { return chrT->overE.sn();   }
                        case AlignBase:   { return chrT->overB.m.sn(); }
                        case AlignIntron: { return chrT->overI.sn();   }
                    }
                }
                
                // Sensitivity at the gene level
                inline double sn(const GeneID &id) const
                {
                    return chrT->geneE.at(id).sn();
                }
                
                // Overall precision
                inline double pc(enum AlignMetrics m) const
                {
                    switch (m)
                    {
                        case AlignExon:   { return chrT->overE.precise(); }
                        case AlignBase:   { return chrT->overB.m.ac();    }
                        case AlignIntron: { return chrT->overI.precise(); }
                    }
                }
            };

            // Analyze a single sample
            static Stats analyze(const FileName &, const Options &o = Options());
        
            // Analyze a single sample
            static Stats analyze(const std::vector<Alignment> &, const Options &o = Options());

            // Analyze multiple replicates
            static std::vector<Stats> analyze(const std::vector<FileName> &, const Options &o = Options());

            // Analyze multiple replicates
            static std::vector<Stats> analyze(const std::vector<std::vector<Alignment>> &,
                                              const Options &o = Options());

            static void report(const FileName &, const Options &o = Options());
            static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif