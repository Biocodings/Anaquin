#ifndef GI_T_ABUND_HPP
#define GI_T_ABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TExpress : public Analyzer
    {
        enum RNALevel
        {
            Gene,
            Isoform
        };

        struct Stats : public LinearStats, public MappingStats
        {
            // The keys depend on whether it's a gene or isoform analysis
            std::map<std::string, Counts> h;
        };

        struct Options : public AnalyzerOptions
        {
            // This's required by gcc...
            Options() {}

            RNALevel level = Isoform;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif
