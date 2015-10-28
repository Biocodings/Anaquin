#ifndef F_NORMAL_HPP
#define F_NORMAL_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FNormal
    {
        typedef FuzzyOptions Options;

        struct Stats : public LinearStats, public FusionStats
        {
            //SequinHist h;
            //Limit ss;
        };

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif