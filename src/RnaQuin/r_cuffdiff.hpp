/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_CUFFDIFF_HPP
#define R_CUFFDIFF_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RCuffdiff
    {
        struct Stats
        {
            struct Data
            {
                ChrID cID;
                
                GeneID gID;
                
                // Not always available
                IsoformID iID;
                
                Fold logF;
                
                // P-value probability
                Probability p;
                
                // Q-value probability
                Probability q;
            };
            
            std::vector<Data> data;
        };
        
        typedef AnalyzerOptions Options;

        static Stats stats (const FileName &, const Options &o = Options());
        static void  report(const FileName &, const Options &o = Options());
    };
}

#endif