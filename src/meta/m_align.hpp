#ifndef M_ALIGN_HPP
#define M_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAlign
    {
        typedef AnalyzerOptions Options;

        struct Stats : public MappingStats
        {
            Sensitivity ss;
            
            // Distribution of the sequins
            SequinHist h = Standard::instance().r_meta.hist();
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif