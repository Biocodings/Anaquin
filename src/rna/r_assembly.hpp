#ifndef GI_R_ASSEMBLY_HPP
#define GI_R_ASSEMBLY_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct RAssembly : RAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats
        {
            Counter cb = RAnalyzer::geneCounter();
            Counter ce = RAnalyzer::isoformCounter();
            Counter ci = RAnalyzer::isoformCounter();
            Counter ct = RAnalyzer::isoformCounter();

            // Performance for each level
            Performance p, pb, pe, pt, pi;
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif