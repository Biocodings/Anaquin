#ifndef R_KALLISTO_HPP
#define R_KALLISTO_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RKallisto
    {
        typedef AnalyzerOptions Options;

        static void analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif