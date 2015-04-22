#ifndef GI_STRUCTURAL_HPP
#define GI_STRUCTURAL_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct StructuralStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct Structural
    {
        enum StructuralLevel
        {
            Base,
            Homozygous,
            Heterzygous
        };

        struct Options : public AnalyzerOptions<StructuralLevel>
        {
            // Empty Implementation
        };

        static StructuralStats analyze(const std::string &file, const Structural::Options &options = Structural::Options());
    };    
}

#endif