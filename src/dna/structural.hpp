#ifndef GI_STRUCTURAL_HPP
#define GI_STRUCTURAL_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct StructuralStats
    {
        // Overall performance
        Confusion m;

        // Overall sensitivity
        Sensitivity s;

        // Base-level performance
        Confusion mb;

        // Base-level sensitivity
        Sensitivity sb;

        // Percentage of variants detected
        double covered;

        // Performance for each genotype
        std::map<Genotype, Confusion> m_gts;

        // Allelle frequencies
        std::map<Sequence, std::pair<Counts, Counts>> afs;
    };

    struct Structural
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static StructuralStats analyze(const std::string &file, const Options &options = Options());
    };    
}

#endif