/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_REPORT_HPP
#define T_REPORT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TReport
    {
        typedef ReportOptions Options;
        
        static void generate(const FileName &file1, const FileName &file2, const Options &o = Options())
        {
            Script::report("TransQuin", file1, file2, o);
        }
    };
}

#endif