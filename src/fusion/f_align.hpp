#ifndef GI_F_ALIGN_HPP
#define GI_F_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FAlign
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : ModelStats
        {
            // Overall performance
            Performance p;

            SequinHist c = Analyzer::histogram(Standard::instance().f_seqs_A);
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif