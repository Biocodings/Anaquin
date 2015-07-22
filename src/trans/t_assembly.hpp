#ifndef GI_T_ASSEMBLY_HPP
#define GI_T_ASSEMBLY_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TAssembly : TAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats
        {
            Counter cb = TAnalyzer::geneCounter();
            Counter ce = TAnalyzer::sequinCounter();
            Counter ci = TAnalyzer::sequinCounter();
            Counter ct = TAnalyzer::sequinCounter();

            // Performance for each level
            Performance p, pb, pe, pt, pi;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif