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
    struct UnknownAlignment
    {
        UnknownAlignment(const std::string &id, const Locus &l) : l(l) {}
        
        // Eg: HISEQ:132:C7F8BANXX:7:1116:11878:9591
        std::string id;
        
        // The position of the alignment
        Locus l;
    };

    class RAlign : public Analyzer
    {
        public:

            typedef AnalyzerOptions Options;

            struct Stats : public AlignmentStats
            {
                struct Data
                {
                    struct AlignLevel
                    {
                        // Number of spliced alignments
                        Counts spliced = 0;
                        
                        // Number of normal alignments
                        Counts normal = 0;
                        
                        Confusion m;
                    };
                    
                    struct IntronLevel
                    {
                        // Unique introns considered FP
                        std::set<Locus> fp;

                        // Confusion for unique introns
                        Confusion m;
                    };
                    
                    struct BaseLevel
                    {
                        // Used for FP aligning outside the reference region
                        std::shared_ptr<MergedInterval> fp;
                    };
                    
                    typedef Confusion ExonLevel;

                    BaseLevel   bLvl;
                    ExonLevel   eLvl;
                    AlignLevel  aLvl;
                    IntronLevel iLvl;

                    std::map<ExonID, Counts> e2r;
                    std::map<IntronID, Counts> i2r;
                    
                    std::map<GeneID, Counts> g2r;
                    std::map<GeneID, Confusion> g2i;
                    
                    // Alignments that have no mapping
                    //std::vector<UnknownAlignment> unknowns;
                };

                std::map<ChrID, Data> data;
                
                MC2Intervals eInters;
                MC2Intervals iInters;

                /*
                 * Synthetic statistics
                 */

                Counts sn = 0;
                Counts ss = 0;
                Confusion sbm, sam, sim, sem;
                
                /*
                 * Genomic statistics
                 */

                Counts gn = 0;
                Counts gs = 0;
                Confusion gbm, gam, gim, gem;
            };

            static Stats analyze(const FileName &, const Options &o = Options());
            static void  report (const FileName &, const Options &o = Options());
    };
}

#endif