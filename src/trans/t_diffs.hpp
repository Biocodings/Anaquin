#ifndef T_DIFFS_HPP
#define T_DIFFS_HPP

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
            Sensitivity ss;

            // The keys depend on whether it's a gene or isoform analysis
            std::map<std::string, Counts> h;
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif