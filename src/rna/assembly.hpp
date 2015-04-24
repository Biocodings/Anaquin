#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "classify.hpp"
#include "r_analyzer.hpp"

namespace Spike
{
    struct AssemblyStats : public AnalyzerStats
    {
        Confusion me;
        Confusion mt;
        Confusion mi;

        Sensitivity se;
        Sensitivity st;
        Sensitivity si;
    };

    struct Assembly : RAnalyzer
    {
        enum Level
        {
            Base,
            Exon,
            Intron,
            Transcripts,
        };

        struct Options : public SingleMixtureOptions<Assembly::Level>
        {
            // Empty Implementation
        };

        static AssemblyStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif