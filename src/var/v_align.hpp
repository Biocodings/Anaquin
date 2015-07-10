#ifndef GI_D_ALIGN_HPP
#define GI_D_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Spike
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
            Counts n_samps = 0;
            
            // Fraction of sequin spiked
            inline Percentage dilution() const { return n_chrT / (n_chrT + n_samps); }

            Performance p;
            
            // Counts for each sequin
            Counter c = DAnalyzer::counterSequins();
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif