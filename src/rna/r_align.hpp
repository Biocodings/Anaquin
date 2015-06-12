#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct RAlign : public RAnalyzer
    {
        struct Stats
        {
            // Total mapped to the in-silico chromosome
            Counts n_chromo = 0;
            
            // Total mapped to the samples
            Counts n_samps = 0;
            
            inline Percentage dilution() const { return n_chromo / (n_chromo + n_samps); }

            // Number of times that each exon is positively identified
            LocusCounter ec = RAnalyzer::exonCounter();
            
            // Number of times that each intron is positively identified
            LocusCounter ic = RAnalyzer::intronCounter();

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
    };
}

#endif