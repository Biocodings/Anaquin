/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef V_KEXPRESS_HPP
#define V_KEXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKExpress
    {
        typedef IndexOptions Options;
        
        struct Stats : public MappingStats, public SequinStats, public LinearStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif