#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "sequins.hpp"
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
        static AssemblyStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
    };
}

#endif