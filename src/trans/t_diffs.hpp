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

        struct Stats : public LinearStats, public MappingStats
        {
            // The keys depend on whether it's a gene or isoform analysis
            std::map<std::string, Counts> h;
        };

        static Stats report(const std::string &, const Options &options = Options());
    };
}

#endif