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
            // Performance at the base level
            Performance bp;
            
            // Performance at the sequin level
            Performance sp;
            
            // Performance for all species at the base level
            std::map<GenomeID, Confusion> base;

            // Reference community
            Intervals<> inters;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif