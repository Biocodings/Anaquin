#ifndef V_VSCAN_HPP
#define V_VSCAN_HPP

#include "data/types.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VVScan
    {
        typedef MappingStats Stats;
        typedef AnalyzerOptions Options;
        
        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif