/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_KDIFF_HPP
#define T_KDIFF_HPP

#include "data/dtest.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TKDiff : public Analyzer
    {
        enum class Software
        {
            Sleuth,
        };
        
        struct Options : public DoubleMixtureOptions
        {
            Options() {}

            // Only Kallisto is supported
            Software soft;
        };

        struct Stats : public MappingStats, public SequinStats, public LinearStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const Options &o);
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
