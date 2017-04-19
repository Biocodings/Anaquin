#ifndef V_ALIGN_HPP
#define V_ALIGN_HPP

#include "data/data.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef AnalyzerOptions Options;
        
        struct Performance
        {
            struct Data
            {
                struct BaseLevel
                {
                    // Used for FP aligning outside the reference region
                    std::shared_ptr<MergedInterval> fp;
                };
                
                struct AlignLevel
                {
                    Confusion m;
                    
                    // Regions to reads (sequins for synthetic)
                    std::map<RegionID, Coverage> r2r;
                };
                
                // Base level for a chromosome
                BaseLevel bLvl;
                
                // Alignment level for a chromosome
                AlignLevel aLvl;
                
                // Overall TP and FP (for each chromosome)
                Counts tp, fp;

                // FP alignments (overlaps)
                std::vector<ReadName> afp;

                std::map<GeneID, Base> lGaps;
                std::map<GeneID, Base> rGaps;
                std::map<GeneID, Base> align;
            };

            std::map<ChrID, Data> data;

            std::map<ChrID, MergedIntervals<>> inters;
            
            /*
             * ----------------- Aggregated statistics -----------------
             */
            
            /*
             * Overall performance at the alignment level
             */
            
            Confusion align;
            
            /*
             * Overall performance at the base level
             */

            Confusion base;
            
            // Sequins to covered
            std::map<SequinID, Base> covered;
            
            // Sequins to length
            std::map<SequinID, Base> length;
            
            // Precision for each region
            std::map<SequinID, Proportion> r2p;
            
            // Sensitivity for each region
            std::map<SequinID, Proportion> r2s;
            
            // Number of mapped alignments
            Counts nMap;
            
            // Number of unmapped alignments
            Counts nNA;
        };
        
        struct Stats
        {
            std::shared_ptr<Performance> endo;
            std::shared_ptr<Performance> seqs;
            
            // Proportion of alignments for endogenous
            inline Proportion pEndo() const
            {
                return (Proportion) endo->nMap / (endo->nMap + seqs->nMap);
            }
            
            // Proportion of alignments for endogenous
            inline Proportion pSeqs() const
            {
                return 1.0 - pEndo();
            }
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o);
        
        static void report(const FileName &, const FileName &, const Options &o = Options());
        
        static void writeSummary(const FileName &,
                                 const FileName &,
                                 const FileName &,
                                 const VAlign::Stats &,
                                 const VAlign::Options &);
        
        static void writeQuins(const FileName &,
                               const VAlign::Stats &,
                               const VAlign::Options &);
        
        static void writeQueries(const FileName &,
                                 const VAlign::Stats &,
                                 const VAlign::Options &);

        static void writeBQuins(const FileName &,
                                const VAlign::Stats &,
                                const VAlign::Options &);
    };
}

#endif
