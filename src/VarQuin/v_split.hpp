#ifndef V_SPLIT_HPP
#define V_SPLIT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSplit
    {
        typedef AnalyzerOptions Options;

        static void report(const FileName &, const Options &);
    };
}

#endif