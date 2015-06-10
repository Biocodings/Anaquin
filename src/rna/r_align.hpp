#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct RAlign : public RAnalyzer
    {
        struct Stats : public AnalyzerStats
        {
            // Number of times that each exon is positively identified
            LocusCounter ec = RAnalyzer::exonCounter();
            
            // Number of times that each intron is positively identified
            LocusCounter ic = RAnalyzer::intronCounter();

            Counter cb = RAnalyzer::geneCounter();
            Counter ce = RAnalyzer::geneCounter();
            Counter ci = RAnalyzer::geneCounter();

            Confusion   mb, me, mi;
            Sensitivity sb, se, si;
        };

        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif