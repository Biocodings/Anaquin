#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RAlignStats : public AnalyzerStats
    {
        Confusion mb;
        Confusion me;
        Confusion mj;

        Sensitivity sb;
        Sensitivity se;
        Sensitivity sj;
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