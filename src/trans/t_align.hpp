#ifndef GI_T_ALIGN_HPP
#define GI_T_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class TAlign : public Analyzer
    {
        public:
            typedef FuzzyOptions Options;

            struct Stats : public MappingStats
            {
                Counts unmapped;
                
                TransRef::GeneHist hb = Standard::instance().r_trans.histGene();
                TransRef::GeneHist he = Standard::instance().r_trans.histGene();
                TransRef::GeneHist hi = Standard::instance().r_trans.histGene();

                // Metrics at various levels
                Performance pb, pe, pi;
            };

            static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif