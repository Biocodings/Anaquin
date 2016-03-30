#ifndef V_EXPRESS_HPP
#define V_EXPRESS_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VExpress
    {
        enum Software
        {
            Kallisto,
        };
        
        struct Options : public AnalyzerOptions
        {
            Software soft;
        };

        struct Stats : public MappingStats, public SequinStats, public LinearStats
        {
            // Absolute detection limit
            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif