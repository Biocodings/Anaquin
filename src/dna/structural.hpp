#ifndef GI_STRUCTURAL_HPP
#define GI_STRUCTURAL_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct StructuralStats : public AnalyzerStats
    {
        Confusion m;
        Percentage covered;
        Percentage efficiency;
    };

    struct Structural
    {
        enum StructuralLevel
        {
            LevelBase,
            LevelHomozygous,
            LevelHeterzygous
        };

        struct Options : public AnalyzerOptions<StructuralLevel>
        {
            // Empty Implementation
        };

        static StructuralStats analyze(const std::string &file, const Structural::Options &options = Structural::Options());
    };    
}

#endif