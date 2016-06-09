#ifndef R_CUFFDIFF_HPP
#define R_CUFFDIFF_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RCuffdiff
    {
        typedef AnalyzerOptions Options;
        typedef MappingStats Stats;

        static Stats stats(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif