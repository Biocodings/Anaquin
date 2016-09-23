/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef V_REPORT_HPP
#define V_REPORT_HPP

#include "tools/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VReport
    {
        typedef ReportOptions Options;

        static void generate(const FileName &file1, const FileName &file2, const Options &o = Options())
        {
//            Script::report("VarQuin", file1, file2, o);
        }
    };
}

#endif