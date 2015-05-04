#ifndef GI_R_ASSEMBLY_HPP
#define GI_R_ASSEMBLY_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RAssemblyStats : public AnalyzerStats
    {
        // Number of times that each exon is positively identified
        LocusCounter e_lc = RAnalyzer::exonCounter();
        
        // Number of times that each intron is positively identified
        LocusCounter i_lc = RAnalyzer::intronCounter();

        // Number of times that each isoform (transcript) is positively identified
        IsoformCounter t_lc = RAnalyzer::isoformCounter();

        Counter cb = RAnalyzer::isoformCounter();
        Counter ce = RAnalyzer::isoformCounter();
        Counter ci = RAnalyzer::isoformCounter();
        Counter ct = RAnalyzer::isoformCounter();

        Confusion mb, me, mt, mi;
        Sensitivity sb, se, st, si;
    };

    struct RAssembly : RAnalyzer
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static RAssemblyStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif