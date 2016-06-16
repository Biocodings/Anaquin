#ifndef R_SLEUTH_HPP
#define R_SLEUTH_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RSleuth
    {
        typedef AnalyzerOptions Options;

        static void analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif