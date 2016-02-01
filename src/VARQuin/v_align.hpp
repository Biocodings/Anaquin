#ifndef V_ALIGN_HPP
#define V_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public AlignmentStats
        {
            Performance p;

            // Distribution for the sequin pairs (reference + variant)
            VarRef::GenoHist h = Standard::instance().r_var.genoHist();
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif