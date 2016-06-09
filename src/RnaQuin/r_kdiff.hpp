/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_KDIFF_HPP
#define R_KDIFF_HPP

#include "stats/analyzer.hpp"
#include "RnaQuin/r_diff.hpp"

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