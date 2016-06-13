#ifndef R_DESEQ2_HPP
#define R_DESEQ2_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RDESeq2
    {
        typedef AnalyzerOptions Options;

        static void analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif