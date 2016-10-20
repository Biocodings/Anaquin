/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef M_DIFF_HPP
#define M_DIFF_HPP

#include "stats/analyzer.hpp"
#include "MetaQuin/m_abund.hpp"

namespace Anaquin
{
    struct MDiff
    {
        struct Stats : public SequinStats, public MappingStats
        {
            MAbund::Stats stats1, stats2;
        };
        
        typedef MAbund::Format Format;
        
        struct Options : public DoubleMixtureOptions
        {
            Format format;
        };

        static Scripts generateRLinear(const FileName &, const Stats &, const Options &);
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
