#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "analyzer.hpp"
#include "sensitivity.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    struct AssemblyStats : public AnalyzerStats
    {
        ConfusionMatrix m_base;
        ConfusionMatrix m_exon;
        ConfusionMatrix m_trans;
        ConfusionMatrix m_intron;
        
        Sensitivity sens_base;
        Sensitivity sens_exon;
        Sensitivity sens_trans;
        Sensitivity sens_intron;
    };

    struct Assembly
    {
        struct AssemblyOptions : public AnalyzerOptions
        {
            // Empty Implementation
        };
        
        static AssemblyStats analyze(const std::string &file, const AssemblyOptions &options = AssemblyOptions());
    };
}

#endif