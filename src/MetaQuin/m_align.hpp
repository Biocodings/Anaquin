#ifndef M_ALIGN_HPP
#define M_ALIGN_HPP

#include "data/data.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAlign
    {
        typedef AnalyzerOptions Options;
        
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
                    std::map<SequinID, Coverage> r2r;
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

        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<Stats> stats;
            
            for (const auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }
            
            return stats;
        }

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
