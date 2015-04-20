#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "classify.hpp"
#include "analyzer.hpp"

namespace Spike
{
    struct AssemblyStats : public AnalyzerStats
    {
        Confusion m_exon;
        Confusion m_trans;
        Confusion m_intron;

        Sensitivity s_exon;
        Sensitivity s_trans;
        Sensitivity s_intron;
    };

    struct Assembly
    {
        enum Mode
        {
            Base,
            Exon,
            Intron,
            Transcripts,
        };

        struct Options : public AnalyzerOptions<Assembly::Mode>
        {
            // Empty Implementation
        };

        static AssemblyStats analyze(const std::string &file, const Assembly::Options &options = Assembly::Options());
    };
}

#endif