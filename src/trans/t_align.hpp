#ifndef T_ALIGN_HPP
#define T_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class TAlign : public Analyzer
    {
        public:
            typedef FuzzyOptions Options;

            struct Stats : public AlignmentStats
            {
                Counts unmapped = 0;

                // Metrics at various levels
                Performance pb, pe, pi;
                
                TransRef::GeneHist hb = Standard::instance().r_trans.histGene();
                TransRef::GeneHist he = Standard::instance().r_trans.histGene();
                TransRef::GeneHist hi = Standard::instance().r_trans.histGene();
            };

            static Stats report(const std::string &, const Options &options = Options());
    };
}

#endif