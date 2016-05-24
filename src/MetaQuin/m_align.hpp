#ifndef M_ALIGN_HPP
#define M_ALIGN_HPP

#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MAlign
    {
        struct Stats : public LinearStats, public AlignmentStats, public SequinStats
        {
            // Empty Implementation
        };

        typedef AnalyzerOptions Options;

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif