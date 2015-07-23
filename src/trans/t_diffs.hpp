#ifndef GI_T_DIFFS_HPP
#define GI_T_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TDiffs : public TAnalyzer
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

        struct Stats : public ModelStats
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &f, const Options &options = Options());
    };
}

#endif