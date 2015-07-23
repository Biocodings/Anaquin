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

        struct Stats : public ModelStats
        {
            // Total mapped to the in-silico chromosome
            Counts n_chrT = 0;
            
            // Total mapped to the samples
            Counts n_genome = 0;

            // Fraction of sequin spiked
            inline Percentage dilution() const { return n_chrT / (n_chrT + n_genome); }

            Performance p;

            // Counts for each sequin
            //SequinHist c = Analyzer::histogram(Standard::instance().v_seqs_A);
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif