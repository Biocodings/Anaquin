#ifndef V_PROCESS_HPP
#define V_PROCESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VProcess
    {
        typedef AnalyzerOptions Options;

        static void report(const FileName &, const Options &);
    };
}

#endif
