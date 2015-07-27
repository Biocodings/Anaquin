#ifndef GI_F_DISCOVER_HPP
#define GI_F_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FExpress
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : ModelStats
        {
            // Overall performance
            Performance p;

            SequinHist h = Analyzer::histogram(Standard::instance().f_seqs_A);
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif