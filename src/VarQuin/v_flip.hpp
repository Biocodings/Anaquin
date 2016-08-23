/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef V_FLIP_HPP
#define V_FLIP_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VFlip
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public MappingStats
        {
            // Low mapping quality
            Counts countLow = 0;
        };

        static Stats analyze(const FileName &,
                             const FileName &,
                             const FileName &,
                             const Options &o = Options());

        static void report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif