#ifndef L_ABUND_HPP
#define L_ABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LAbund
    {
        typedef AnalyzerOptions Options;

        struct Stats : public LinearStats, public AlignmentStats
        {
            // Histogram expected
            std::map<SequinID, Coverage> expect;

            // Histogram before normalization and adjustment
            std::map<SequinID, Counts> measured;

            // Histogram after normalization but before adjustment
            std::map<SequinID, Coverage> normalized;

            // Histogram after adjustment
            std::map<SequinID, Coverage> adjusted;
            
            // Adjusted abundance for each sequin
            std::map<BaseID, Coverage> sequinAdjusted;

            // Expected size of the library
            Counts expTotal = 0;

            // Measured size of the library
            Counts obsTotal = 0;

            SequinHist h = Standard::instance().r_lad.hist();
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif