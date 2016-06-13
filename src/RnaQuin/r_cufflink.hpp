#ifndef R_CUFFLINK_HPP
#define R_CUFFLINK_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RCufflink
    {
        typedef AnalyzerOptions Options;

        static void analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif