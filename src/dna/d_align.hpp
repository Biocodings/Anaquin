#ifndef GI_D_ALIGN_HPP
#define GI_D_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct DAlign
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };
        
        struct Stats : public ModelStats
        {
            // Total mapped to the in-silico chromosome
            Counts n_chromo = 0;
            
            // Total mapped to the samples
            Counts n_samps = 0;
            
            inline Percentage dilution() const { return n_chromo / (n_chromo + n_samps); }
            
            Performance p;
            Counter c = DAnalyzer::counterSequins();
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif