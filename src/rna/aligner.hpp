#ifndef GI_ALIGNER_HPP
#define GI_ALIGNER_HPP

#include "analyzer.hpp"
#include "sensitivity.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    struct AlignerStats : public AnalyzerStats
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
        
        struct AlignerOptions : public AnalyzerOptions
        {
            AlignerMode mode;
        };

        static AlignerStats analyze(const std::string &file, const AlignerOptions &options = AlignerOptions());
    };
}

#endif