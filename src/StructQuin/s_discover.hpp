#ifndef S_DISCOVER_HPP
#define S_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct SDiscover : public Analyzer
    {
        struct Stats : public MappingStats
        {

        };
        
        typedef AnalyzerOptions Options;

        static Stats analyze(const FileName &, const Options &);
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
