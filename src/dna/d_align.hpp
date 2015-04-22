#ifndef GI_D_ALIGN_HPP
#define GI_D_ALIGN_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct DAlignStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct DAlign
    {
        enum Mode
        {
            Base,
        };

        struct Options : public AnalyzerOptions<DAlign::Mode>
        {
            // Empty Implementation
        };

        static DAlignStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif