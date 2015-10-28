#ifndef F_DIFF_HPP
#define F_DIFF_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FDiff
    {
        typedef FuzzyOptions Options;

        struct Stats : public LinearStats, public FusionStats
        {
            Limit ss;

            // Sequin distribution
            SequinHist h = Standard::instance().r_fus.hist();;
        };

        static Stats stats(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif