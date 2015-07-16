#ifndef GI_R_ABUND_HPP
#define GI_R_ABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RAbund : public RAnalyzer
    {
        struct Stats : public ModelStats
        {
            SequinCounter c = RAnalyzer::sequinCounter();
        };

        struct Options : public SingleMixtureOptions
        {
            //Options() {}
            RNALevel level = Isoform;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif