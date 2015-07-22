#ifndef GI_T_ALIGN_HPP
#define GI_T_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class TAlign : public TAnalyzer
    {
        public:
            struct Stats
            {
                // Total mapped to the in-silico chromosome
                Counts n_chrT = 0;

                // Total mapped to the human genome
                Counts n_genome = 0;

                // Tracker for each exon for each gene
                GeneTracker g_exon_tracker = TAnalyzer::geneTracker();

                // Tracker for each intron for each gene
                GeneTracker g_intron_tracker = TAnalyzer::geneTracker();

                // Fraction of sequin spiked
                inline Percentage dilution() const { return n_chrT / (n_chrT + n_genome); }

                Counter cb = TAnalyzer::geneCounter();
                Counter ce = TAnalyzer::geneCounter();
                Counter ci = TAnalyzer::geneCounter();

                Performance pb, pe, pi;
            };

            struct Options : public SingleMixtureOptions
            {
                // Empty Implementation
            };

            static Stats analyze(const std::string &, const Options &options = Options());

        private:
            void reportGeneral(const std::string &, const Stats &, const Options &);
    };
}

#endif