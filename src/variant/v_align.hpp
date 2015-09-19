#ifndef GI_V_ALIGN_HPP
#define GI_V_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef SingleMixtureOptions Options;
        
        struct Stats : public AlignmentStats
        {
            Performance p;

            // Distribution for the sequin pairs (reference + variant)
            VarRef::PairHist h = Standard::instance().r_var.pairHist();
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif