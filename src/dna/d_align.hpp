#ifndef GI_D_ALIGN_HPP
#define GI_D_ALIGN_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct DAlignStats
    {
        // Empty Implementation
    };

    struct DAlign
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static DAlignStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif