#ifndef V_REPORT_HPP
#define V_REPORT_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VReport
    {
        typedef ReportOptions Options;

        static void generate(const FileName &file1, const FileName &file2, const Options &o = Options())
        {
            Script::report("VarQuin", file1, file2, o);
        }
    };
}

#endif