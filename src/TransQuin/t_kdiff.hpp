/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_KDIFF_HPP
#define T_KDIFF_HPP

#include "stats/analyzer.hpp"
#include "TransQuin/t_diff.hpp"

namespace Anaquin
{
    struct TKDiff : public Analyzer
    {
        typedef IndexOptions Options;
        typedef TDiff::Stats Stats;
        
        static Stats analyze(const std::vector<FileName> &,
                             const std::vector<FileName> &,
                             const std::vector<FileName> &,
                             const std::vector<FileName> &,
                             const Options &o);

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif