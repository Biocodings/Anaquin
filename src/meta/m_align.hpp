#ifndef M_ALIGN_HPP
#define M_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAlign
    {
        typedef AnalyzerOptions Options;

        struct Stats : public AlignmentStats
        {
            Performance p;

            // Distribution of the sequins
            SequinHist h = Standard::instance().r_meta.hist();
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif