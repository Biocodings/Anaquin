/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_COVERAGE_HPP
#define R_COVERAGE_HPP

#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TCoverage
    {
        typedef AnalyzerOptions Options;        
        typedef CoverageTool::Stats Stats;

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif