#ifndef M_ABUND_HPP
#define M_ABUND_HPP

#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MAbund
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