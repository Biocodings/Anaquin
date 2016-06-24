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
            };

            std::map<ChrID, Data> data;

            // Histogram for all chromosomes
            std::map<ChrID, Hist> hist;
            
            std::map<ChrID, MergedIntervals<>> inters;
            
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