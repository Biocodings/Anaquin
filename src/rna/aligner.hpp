#ifndef GI_ALIGNER_HPP
#define GI_ALIGNER_HPP

#include "sensitivity.hpp"
#include "parsers/parser.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    struct AlignerStats
    {
        // Percentage of reads aligned with the reference chromosome
        Percentage pr;
        
        // Percentage of reads aligned with the query samples
        Percentage pq;
        
        // Total number of reads aligned
        Reads n = 0;
        
        // Number of reads aligned to the chromosome
        Reads nr = 0;
        
        // Number of reads aligned to the real sample
        Reads nq = 0;

        // The lowest detectable abundance in the experiment
        Sensitivity sens;

        // Amount of sequins detected
        Percentage dilution;
        
        ConfusionMatrix m;        
    };

    struct Aligner
    {
        enum AlignerMode
        {
            BaseAlign,
            ExonAlign,
            SpliceAlign,
        };
        
        struct AlignerOptions : public ParserOptions
        {
            AlignerMode mode;
        };

        static AlignerStats analyze(const std::string &file, const AlignerOptions &options = AlignerOptions());
    };
}

#endif