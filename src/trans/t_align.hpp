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
                Counts unmapped = 0;

                // Metrics at various levels
                Performance pb, pe, pi;

                // Intervals for exons in TransQuin reference
                Intervals<TransRef::ExonInterval> exonInters;
                
                // Intervals for introns in TransQuin reference
                Intervals<TransRef::IntronInterval> intronInters;
            };

            static Stats report(const std::string &, const Options &options = Options());
    };
}

#endif