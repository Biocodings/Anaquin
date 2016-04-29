#ifndef F_PARENT_HPP
#define F_PARENT_HPP

#include "stats/analyzer.hpp"
#include "FusQuin/FUSQuin.hpp"

namespace Anaquin
{
    struct FParent
    {
        struct Options : public FuzzyOptions
        {
            FusionCaller caller;
        };

        struct Stats : public FusionStats, public LinearStats, public SequinStats
        {
            // Absolute detection limit
            Limit absolute;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif