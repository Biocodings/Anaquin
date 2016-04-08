/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_KEXPRESS_HPP
#define T_KEXPRESS_HPP

#include "stats/analyzer.hpp"
#include "TransQuin/t_express.hpp"

namespace Anaquin
{
    struct TKExpress : public Analyzer
    {
        typedef IndexOptions Options;
        typedef TExpress::Stats Stats;

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
