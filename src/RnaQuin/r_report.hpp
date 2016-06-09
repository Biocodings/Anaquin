/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_REPORT_HPP
#define R_REPORT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TReport
    {
        typedef ReportOptions Options;
        
        /*
         * Generate a report for paired-end sequence files. This'll only work for expression analysis.
         */
        
        static void generate(const FileName &file1, const FileName &file2, const Options &o = Options())
        {
            Script::report("TransQuin", file1, file2, o);
        }
        
        /*
         * Generate a report for metadata. This'll work for both expression and differential analysis.
         */
        
        static void generate(const FileName &meta, const Options &o = Options())
        {
            Script::report("TransQuin", meta, meta, o);
        }
    };
}

#endif