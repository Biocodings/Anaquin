#ifndef GI_T_ALIGN_HPP
#define GI_T_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class TAlign : public Analyzer
    {
        public:
            typedef FuzzyOptions Options;

            struct Stats
            {
                // Total mapped to the in-silico chromosome
                Counts n_chrT = 0;

                // Total mapped to the human genome
                Counts n_genome = 0;

                // Fraction of sequin spiked
                inline Percentage dilution() const { return n_chrT / (n_chrT + n_genome); }

                BaseHist hb = Analyzer::baseHist();
                BaseHist he = Analyzer::baseHist();
                BaseHist hi = Analyzer::baseHist();

                Performance pb, pe, pi;
            };

            static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif