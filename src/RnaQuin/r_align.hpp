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
                    
                    struct BaseLevel
                    {
                        Base fp = 0;
                    };
                    
                    typedef Confusion ExonLevel;
                    typedef Confusion IntronLevel;

                    BaseLevel   bLvl;
                    ExonLevel   eLvl;
                    AlignLevel  aLvl;
                    IntronLevel iLvl;

                    std::map<ExonID, Counts> e2r;
                    std::map<IntronID, Counts> i2r;

                    // Alignments that have no mapping
                    //std::vector<UnknownAlignment> unknowns;
                };

                std::map<ChrID, Data> data;
                
                MC2Intervals eInters;
                MC2Intervals iInters;

                // Synthetic Statistics
                Confusion sbm, sam, sim, sem;
                
                // Genomic Statistics
                Confusion gbm, gam, gim, gem;
            };

            static Stats analyze(const FileName &, const Options &o = Options());
            static void  report (const FileName &, const Options &o = Options());
    };
}

#endif