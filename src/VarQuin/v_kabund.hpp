#ifndef V_KABUND_HPP
#define V_KABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKAbund
    {
        enum class Software
        {
            Salmon
        };

        enum class Mode
        {
            CNVLad,
            AFLad,
            ConLad,
        };
        
        struct Options : public AnalyzerOptions
        {
            Mode mode;
            Software soft = Software::Salmon;
        };
        
        struct Stats : public LimitStats, public SequinStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const Options &o);

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
