#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "analyzer.hpp"
#include "confusion.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct AssemblyStats : public AnalyzerStats
    {
        Confusion m_base;
        Confusion m_exon;
        Confusion m_trans;
        Confusion m_intron;
        
        Sensitivity sens_base;
        Sensitivity sens_exon;
        Sensitivity sens_trans;
        Sensitivity sens_intron;
    };

    struct Assembly
    {
        enum AssemblyMode
        {
            Assembly_Base,
        };

        struct Options : public AnalyzerOptions<AssemblyMode>
        {
            // Empty Implementation
        };

        static AssemblyStats analyze(const std::string &file, const Assembly::Options &options = Assembly::Options());
    };
}

#endif