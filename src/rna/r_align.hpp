#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "r_analyzer.hpp"

namespace Spike
{
    struct RAlignStats : public AnalyzerStats
    {
        Counter     cb, ce, ci;
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