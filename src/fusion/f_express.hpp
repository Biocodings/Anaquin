#ifndef F_EXPRESS_HPP
#define F_EXPRESS_HPP

#include "fusion/FUSQUin.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FExpress
    {
        struct Options : public FuzzyOptions
        {
            FusionCaller caller;
        };

        struct Stats : public FusionStats, public LinearStats, public SequinStats
        {
            // Absolute detection limit
            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif