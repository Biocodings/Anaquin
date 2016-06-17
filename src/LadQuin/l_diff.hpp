/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef L_DIFF_HPP
#define L_DIFF_HPP

#include "stats/analyzer.hpp"
#include "LadQuin/l_norm.hpp"

namespace Anaquin
{
    struct LDiff
    {
        typedef AnalyzerOptions Options;

        struct Stats : public AlignmentStats
        {
            struct Data : public LinearStats
            {
                // Empty Implementation
            };
            
            std::map<ChrID, Data> data;
            
            LNorm::Stats a, b;
            
            // Sensitivity at the joined level
            Limit s_joined;

            // Sensitivity at the unjoined level
            Limit ss;

            // Histogram at the unjoined level
            SequinHist h = Standard::instance().r_lad.hist();
            
            // Histogram at the joined level
            LadderRef::JoinHist h_joined  = Standard::instance().r_lad.joinHist();
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());        
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif