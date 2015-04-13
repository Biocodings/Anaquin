#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "confusion_matrix.hpp"

namespace Spike
{
    struct AssemblyStats
    {
        ConfusionMatrix base;
        ConfusionMatrix exon;
        ConfusionMatrix intron;
    };
    
    struct Assembly
    {
        static AssemblyStats analyze(const std::string &file);
    };
}

#endif