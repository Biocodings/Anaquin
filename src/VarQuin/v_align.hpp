#ifndef V_ALIGN_HPP
#define V_ALIGN_HPP

#include "data/data.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        struct Options : public AnalyzerOptions
        {
            FileName gBed, sBed;
        };
        
        struct Stats : public AlignmentStats
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
                std::vector<ReadID> afp;

                std::map<GeneID, Base> lGaps;
                std::map<GeneID, Base> rGaps;
                std::map<GeneID, Base> align;
            };

            std::map<ChrID, Data> data;

            std::map<ChrID, MergedIntervals<>> inters;
            
            /*
             * Aggregated statistics
             */
            
            /*
             * Alignment level for synthetic and genome
             */
            
            typedef Confusion AlignLevel;
            
            AlignLevel sa, ga;
            
            /*
             * Base level statistics for synthetic and genome
             */

            Confusion sb;
            Confusion gb;
            
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
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());        
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
