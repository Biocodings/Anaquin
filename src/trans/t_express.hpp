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

        struct Stats : public ModelStats
        {
            SequinHist c = Analyzer::seqHist();
        };

        struct Options : public SingleMixtureOptions
        {
            // This's required by gcc...
            Options() {}

            RNALevel level = Isoform;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif
