#ifndef GI_DENOVO_HPP
#define GI_DENOVO_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct DenovoStats
    {
        // Empty Implementation
    };

    struct Denovo
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static DenovoStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif