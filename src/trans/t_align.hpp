#ifndef T_ALIGN_HPP
#define T_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class TAlign : public Analyzer
    {
        public:
            typedef AnalyzerOptions Options;

            struct Stats : public AlignmentStats
            {
                // Overall performance at various levels
                Performance pb, pe, pi;

                // Individual performance for each sequin
                std::map<GeneID, Confusion> sb, se, si;

                // Intervals for exons in TransQuin
                Intervals<TransRef::ExonInterval> eInters;

                // Intervals for introns in TransQuin
                Intervals<TransRef::IntronInterval> iInters;

                // Sequins that have failed to be detected
                std::map<SequinID, MissingSequin> missings;

                // Alignments that have no mapping
                std::vector<UnknownAlignment> unknowns;
            };

            static Stats stats (const FileName &, const Options &options = Options());
            static void  report(const FileName &, const Options &options = Options());
    };
}

#endif