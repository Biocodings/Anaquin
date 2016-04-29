#ifndef F_FUSION_HPP
#define F_FUSION_HPP

#include "stats/analyzer.hpp"
#include "FusQuin/FUSQuin.hpp"

namespace Anaquin
{
    struct FFusion
    {
        struct Options : public FuzzyOptions
        {
            FusionCaller soft;
        };

        struct Stats : public FusionStats, public LinearStats, public SequinStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif