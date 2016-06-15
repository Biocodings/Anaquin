/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_ALIGN_HPP
#define R_ALIGN_HPP

#include "stats/analyzer.hpp"
#include "data/alignment.hpp"

namespace Anaquin
{
    typedef std::map<BinID, Counts> BinCounts;

    struct UnknownAlignment
    {
        UnknownAlignment(const std::string &id, const Locus &l) : l(l) {}
        
        // Eg: HISEQ:132:C7F8BANXX:7:1116:11878:9591
        std::string id;
        
        // The position of the alignment
        Locus l;
    };

    struct CountProp
    {
        CountProp() {}
        CountProp(Counts i, Counts n) : i(i), n(n) {}
        
        inline operator std::string() const
        {
            return (boost::format("%1% (%2%)") % i % n).str();
        }
        
        inline Proportion percent() const
        {
            assert(i <= n);
            return static_cast<double>(i) / n;
        }
        
        // Relevant number of counts
        Counts i = 0;
        
        // Total number of counts
        Counts n = 0;
    };

    /*
     * Represents something that is missing or undetected. It could be an exon, intron, isoform, gene etc.
     */
    
    struct Missing
    {
        Missing(const GenericID &id) : id(id) {}
        
        inline bool operator==(const Missing &m) const { return id == m.id; }
        inline bool operator< (const Missing &m) const { return id <  m.id; }
        
        const GenericID id;
    };

    class RAlign : public Analyzer
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
                 * Metric statistics (eg: exons and introns)
                 */
                
                Counts lTP = 0;
                Counts lNR = 0;

                inline Counts aNQ() const { return aTP + aFP; }
                inline Counts lFN() const { return lNR - lTP; }

                inline double sn() const
                {
                    return (lTP + lFN()) ? static_cast<double>(lTP) / (lTP + lFN()) : NAN;
                }

                inline double pc() const
                {
                    return (aTP + aFP) ? static_cast<double>(aTP) / (aTP + aFP) : NAN;
                }
            };

            struct Stats : public AlignmentStats
            {
                enum class AlignMetrics
                {
                    AlignExon,
                    AlignIntron,
                    AlignBase
                };

                enum class MissingMetrics
                {
                    MissingExon,
                    MissingIntron,
                    MissingGene
                };

                struct Data
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
                    
                    // Number of spliced alignments
                    Counts n_spliced = 0;
                    
                    // Number of non-spliced alignments
                    Counts n_normal = 0;

                    MergedConfusion overE, overI;
                    
                    // Overall performance at the base level
                    Performance overB;
                    
                    Hist histE, histI;
                    
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

                std::map<ChrID, Data> data;
                
                inline Counts countSpliceSyn() const
                {
                    return count(data, [&](const ChrID &cID, const Data &x)
                    {
                        return Standard::isSynthetic(cID) ? x.n_spliced : 0;
                    });
                }

                inline Counts countSpliceGen() const
                {
                    return count(data, [&](const ChrID &cID, const Data &x)
                    {
                        return !Standard::isSynthetic(cID) ? x.n_spliced : 0;
                    });
                }
                
                inline Counts countNormalSyn() const
                {
                    return count(data, [&](const ChrID &cID, const Data &x)
                    {
                        return Standard::isSynthetic(cID) ? x.n_normal : 0;
                    });
                }
                
                inline Counts countNormalGen() const
                {
                    return count(data, [&](const ChrID &cID, const Data &x)
                    {
                        return !Standard::isSynthetic(cID) ? x.n_normal : 0;
                    });
                }
                
                inline CountProp countMiss(const ChrID &cID, MissingMetrics m) const
                {
                    switch (m)
                    {
                        case MissingMetrics::MissingGene:
                        {
                            return CountProp(data.at(cID).missG.size(), data.at(cID).histE.size());
                        }
                            
                        case MissingMetrics::MissingExon:
                        {
                            return CountProp(data.at(cID).missE.size(), data.at(cID).eContains.size());
                        }

                        case MissingMetrics::MissingIntron:
                        {
                            return CountProp(data.at(cID).missI.size(), data.at(cID).iContains.size());
                        }
                    }
                }

                inline double missProp(const ChrID &cID, MissingMetrics m) const
                {
                    return countMiss(cID, m).percent();
                }
                
                // Overall sensitivity
                inline double sn(const ChrID &cID, AlignMetrics m) const
                {
                    switch (m)
                    {
                        case AlignMetrics::AlignExon:   { return data.at(cID).overE.sn();   }
                        case AlignMetrics::AlignBase:   { return data.at(cID).overB.m.sn(); }
                        case AlignMetrics::AlignIntron: { return data.at(cID).overI.sn();   }
                    }
                }
                
                // Sensitivity at the gene level
                inline Proportion sn(const ChrID &cID, const GeneID &id) const
                {
                    return data.at(cID).geneE.at(id).sn();
                }

                // Overall precision
                inline Proportion pc(const ChrID &cID, AlignMetrics m) const
                {
                    switch (m)
                    {
                        case AlignMetrics::AlignExon:   { return data.at(cID).overE.pc();   }
                        case AlignMetrics::AlignBase:   { return data.at(cID).overB.m.pc(); }
                        case AlignMetrics::AlignIntron: { return data.at(cID).overI.pc();   }
                    }
                }
                
                /*
                 * Synthetic Statistics
                 */
                
                Proportion s_esn, s_epc;
                Proportion s_isn, s_ipc;
                Proportion s_bsn, s_bpc;
                Proportion s_ems, s_ims, s_gms;
                
                // Number of reads mapped to each sequin gene (eg: R1_1)
                std::map<SequinID, Counts> s2r;
                
                /*
                 * Genomic Statistics
                 */
                
                CountProp g_ems, g_ims, g_gms;

                Proportion g_esn = NAN;
                Proportion g_epc = NAN;
                Proportion g_isn = NAN;
                Proportion g_ipc = NAN;
                Proportion g_bsn = NAN;
                Proportion g_bpc = NAN;
            };

            static Stats analyze(const FileName &, const Options &o = Options());
            static Stats analyze(const std::vector<Alignment> &, const Options &o = Options());

            static std::vector<Stats> analyze(const std::vector<std::vector<Alignment>> &aligns, const Options &o = Options())
            {
                std::vector<RAlign::Stats> stats;
                
                for (const auto &align : aligns)
                {
                    stats.push_back(analyze(align, o));
                }
                
                return stats;            
            }

            static void report(const FileName &, const Options &o = Options());
    };
}

#endif