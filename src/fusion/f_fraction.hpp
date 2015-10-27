#ifndef F_FRACTION_HPP
#define F_FRACTION_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FFraction
    {
        typedef FuzzyOptions Options;

        struct Stats : public LinearStats, public FusionStats
        {
            SequinHist h;
            Limit ss;
        };

        static Stats stats(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif