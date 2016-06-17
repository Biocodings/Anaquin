/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef F_DIFF_HPP
#define F_DIFF_HPP

#include "stats/analyzer.hpp"
#include "FusQuin/f_normal.hpp"
#include "FusQuin/f_fusion.hpp"

namespace Anaquin
{
    struct FDiff
    {
        struct Options : public FuzzyOptions
        {
            FusionCaller soft;
        };

        struct Stats : public LinearStats, public FusionStats, public SequinStats
        {
            // Statistics for the normal genes
            FNormal::Stats normal;

            // Statistics for the fusion genes
            FFusion::Stats fusion;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());        
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif