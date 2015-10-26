#ifndef L_DIFFS_HPP
#define L_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LDiffs
    {
        typedef AnalyzerOptions Options;

        struct Stats : public LinearStats, public AlignmentStats
        {
            // Sensitivity at the joined level
            Limit s_joined;
            
            // Histogram at the joined level
            LadderRef::JoinHist h_joined  = Standard::instance().r_lad.joinHist();
            
            // Sensitivity at the unjoined level
            Limit ss;
            
            // Histogram at the unjoined level
            SequinHist h = Standard::instance().r_lad.hist();
        };

        static Stats report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif