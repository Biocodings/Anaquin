#ifndef F_NORMAL_HPP
#define F_NORMAL_HPP

#include "stats/analyzer.hpp"
#include "FusQuin/FUSQuin.hpp"

namespace Anaquin
{
    struct FNormal
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