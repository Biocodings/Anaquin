#ifndef F_DIFF_HPP
#define F_DIFF_HPP

#include "stats/analyzer.hpp"
#include "fusion/f_express.hpp"

namespace Anaquin
{
    struct FDiff
    {
        typedef FuzzyOptions Options;

        struct Stats : public LinearStats, public FusionStats, public SequinStats
        {
            // Statistics for the normal genes
            FExpress::Stats normal;

            // Statistics for the fusion genes
            FExpress::Stats fusion;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());        
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif