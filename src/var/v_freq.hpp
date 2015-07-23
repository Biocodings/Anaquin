#ifndef GI_V_FREQ_HPP
#define GI_V_FREQ_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VFreq
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public ModelStats
        {
            // Overall performance
            Confusion m;

            // The proportion of variations with alignment coverage
            double covered;

            // Measure of variant detection independent to sequencing depth or coverage
            double efficiency;

            LocusHist h = Analyzer::histogram__(Standard::instance().v_vars);
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif