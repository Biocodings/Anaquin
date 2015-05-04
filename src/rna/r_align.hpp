#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RAlignStats : public AnalyzerStats
    {
        /*
         * Counters for various levels and purposes. For example, e_lc is the counter for number of times that each exon
         * positively identified.
         */

        LocusCounter e_lc = RAnalyzer::exonCounter();
        LocusCounter i_lc = RAnalyzer::intronCounter();

        Counter cb = RAnalyzer::geneCounter();
        Counter ce = RAnalyzer::geneCounter();
        Counter ci = RAnalyzer::geneCounter();
        
        Confusion   mb, me, mi;
        Sensitivity sb, se, si;
    };

    struct RAlign : public RAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static RAlignStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif