#ifndef GI_V_ALIGN_HPP
#define GI_V_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public LinearStats
        {
            // Total mapped to the in-silico chromosome
            Counts n_chrT = 0;

            // Total mapped to the samples
            Counts n_genome = 0;

            // Fraction of sequin spiked
            inline Percentage dilution() const { return n_chrT / (n_chrT + n_genome); }

            Performance p;

            // Counts for each sequin
            SequinHist h_seq;//TODO = Analyzer::seqHist();
            
            // Counts for each base
            SequinHist h_base;//TODO = Analyzer::baseHist();
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif