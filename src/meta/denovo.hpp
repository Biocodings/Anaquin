#ifndef GI_DENOVO_HPP
#define GI_DENOVO_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct DenovoStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct Denovo
    {
        enum DenovoLevel
        {
            Base,
        };

        struct Options : public AnalyzerOptions<DenovoLevel>
        {
            // Empty Implementation
        };

        static DenovoStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif