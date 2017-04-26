#ifndef V_KEXPRESS_HPP
#define V_KEXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKExpress
    {
        enum class Software
        {
            Salmon
        };

        struct Options : public AnalyzerOptions
        {
            Software soft;
        };
        
        struct Stats : public LimitStats, public SequinStats
        {
        };

        static Stats analyze(const FileName &, const Options &o);

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
