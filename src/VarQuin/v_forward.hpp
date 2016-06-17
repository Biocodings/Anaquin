/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef V_FORWARD_HPP
#define V_FORWARD_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VForward
    {
        typedef AnalyzerOptions Options;
        
        static void analyze(const FileName &,
                            const FileName &,
                            const FileName &,
                            const Options &o = Options());
        static void report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif