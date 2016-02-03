#ifndef F_NORMAL_HPP
#define F_NORMAL_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FNormal
    {
        typedef FuzzyOptions Options;

        struct Stats
        {
            struct Data : public LinearStats, public FusionStats
            {
                // Detection limit
                Limit ss;
                
                // Sequin distribution
                SequinHist h = Standard::instance().r_fus.hist();
            };

            std::map<ChromoID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif