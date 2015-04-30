#ifndef GI_R_ASSEMBLY_HPP
#define GI_R_ASSEMBLY_HPP

#include "classify.hpp"
#include "r_analyzer.hpp"

namespace Spike
{
    struct RAssemblyStats : public AnalyzerStats
    {
        Confusion mb;
        Confusion me;
        Confusion mt;
        Confusion mi;

        Sensitivity sb;
        Sensitivity se;
        Sensitivity st;
        Sensitivity si;
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