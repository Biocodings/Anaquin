#ifndef V_ALIGN_HPP
#define V_ALIGN_HPP

#include "data/types.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public AlignmentStats
        {
            /*
             * Returns the overall sensitivity, defined by the number of bases aligned
             * with respect to the reference.
             */

            //inline Proportion sn(const ChrID &cID) const
            //{
            //    return static_cast<Proportion>(covered(cID)) / length(cID);
            //}

            // Returns the individual sensitivity
//            inline Proportion sn(const ChrID &cID, const SequinID &sID) const
//            {
//                // How many bases covered?
//                const auto covered = data.at(cID).covered.at(sID);
//                
//                // How long is the reference?
//                const auto length  = data.at(cID).length.at(sID);
//
//                return static_cast<Proportion>(covered) / length;
//            }
//
//            // Returns the overall precision
//            inline Proportion pc(const ChrID &cID) const
//            {
//                const auto &tp = data.at(cID).tp;
//                const auto &fp = data.at(cID).fp;
//
//                return static_cast<Proportion>(tp) / (tp + fp);
//            }

            //inline Base length(const ChrID &cID)  const { return sum(data.at(cID).length);  }
            //inline Base covered(const ChrID &cID) const { return sum(data.at(cID).covered); }
            
            struct Data
            {
                // Overall TP and FP (for each chromosome)
                Counts tp, fp;

                // FP alignments (overlaps)
                std::vector<ReadID> afp;

                std::map<GeneID, Base> lGaps;
                std::map<GeneID, Base> rGaps;
                std::map<GeneID, Base> align;

                // Positions that have no alignment
                std::set<Locus> gaps;
                
//                /*
//                 * Statistics for genes (eg: sequins)
//                 */
//                
//                std::map<GeneID, Counts> gtp;
//                std::map<GeneID, Counts> gfp;
            };

            std::map<ChrID, Data> data;

            // Histogram for all chromosomes
            std::map<ChrID, Hist> hist;
            
            // Intervals for all chromosomes
//            std::map<ChrID, Intervals<>> inters;
            
            std::map<ChrID, MergedIntervals<>> inters;
            
            /*
             * Genomic statistics
             */
            
//            // Total covered by all genes
//            Base gc = 0;
//            
//            // Total length by all genes
//            Base gl = 0;
//            
//            // Number of TP for genomic genes
//            Counts gtp = 0;
//            
//            // Number of FP for genomic genes
//            Counts gfp = 0;
//            
//            // Overall sensitivity for all genomic genes
//            inline Proportion gsn() const
//            {
//                return static_cast<Proportion>(gc) / (gc + gl);
//            }
//            
//            // Overall precision for all genomic genes
//            inline Proportion gpc() const
//            {
//                return static_cast<Proportion>(gtp) / (gtp + gfp);
//            }

            // Sensitivty and precision for the synthetic
            Proportion ssn, spc;
            
            // Sensitivity and precision for the genome
            Proportion gsn, gpc;
            
            // Sequins to covered (synthetic)
            std::map<SequinID, Base> s2c;
            
            // Sequins to length (synthetic)
            std::map<SequinID, Base> s2l;
            
            // Genes to covered (genome)
            std::map<GeneID, Base> g2c;
            
            // Genes to length (genome)
            std::map<GeneID, Base> g2l;
            
            // Genes to precision
            std::map<GeneID, Proportion> g2p;
            
            // Sequins to sensitivity
            std::map<GeneID, Proportion> g2s;

            // Genes to reads
            std::map<SequinID, Coverage> g2r;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif