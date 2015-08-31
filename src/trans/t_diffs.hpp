#ifndef GI_T_DIFFS_HPP
#define GI_T_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TDiffs : public Analyzer
    {
        enum RNALevel
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            RNALevel level;
        };

        struct Stats : public LinearStats
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif