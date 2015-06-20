#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    class RAlign : public RAnalyzer
    {
        public:
            struct Stats
            {
                // Total mapped to the in-silico chromosome
                Counts n_chrT = 0;

                // Total mapped to the samples
                Counts n_samps = 0;

                // Tracker for each exon for each gene
                GeneTracker g_exon_tracker = RAnalyzer::geneTracker();

                // Tracker for each intron for each gene
                GeneTracker g_intron_tracker = RAnalyzer::geneTracker();

                // Fraction of sequin spiked
                inline Percentage dilution() const { return n_chrT / (n_chrT + n_samps); }

                Counter cb = RAnalyzer::geneCounter();
                Counter ce = RAnalyzer::geneCounter();
                Counter ci = RAnalyzer::geneCounter();

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