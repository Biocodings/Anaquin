/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_GENOME_HPP
#define R_GENOME_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class RGenome : public Analyzer
    {
        public:
            typedef AnalyzerOptions Options;
            static void report(const FileName &, const Options &o = Options());
    };
}

#endif