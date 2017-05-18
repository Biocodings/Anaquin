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
	    Options() : mode(Mode::CNVLad), soft(Software::Salmon) {}
            Mode mode;
            Software soft;
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