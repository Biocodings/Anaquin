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
            // Overall performance
            std::map<PerfLevel, Performance> p;

            // Performance for all species at the base level
            std::map<GenomeID, Performance> base;

            // Performance for all species at the sequin level
            std::map<GenomeID, Performance> seq;
            
            // Reference species
            Intervals inters;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif