#ifndef GI_ALIGNER_HPP
#define GI_ALIGNER_HPP

#include "sensitivity.hpp"
#include "parsers/parser.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    struct AlignerStats : public ParserStats
    {
        // The lowest detectable abundance in the experiment
        Sensitivity sens;

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