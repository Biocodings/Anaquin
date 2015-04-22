#ifndef GI_M_ALIGN_HPP
#define GI_M_ALIGN_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct MAlignStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct MAlign
    {
        enum Mode
        {
            Base,
        };

        struct Options : public AnalyzerOptions<MAlign::Mode>
        {
            // Empty Implementation
        };

        static MAlignStats analyze(const std::string &file, const MAlign &options = MAlign());
    };
}

#endif