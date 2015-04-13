#ifndef GI_ASSEMBLY_HPP
#define GI_ASSEMBLY_HPP

#include "parsers/parser.hpp"
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
        struct AssemblyOptions : public ParserOptions
        {
            // Empty Implementation
        };
        
        static AssemblyStats analyze(const std::string &file, const AssemblyOptions &options = AssemblyOptions());
    };
}

#endif