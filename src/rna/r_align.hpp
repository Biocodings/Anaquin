#ifndef GI_R_ALIGN_HPP
#define GI_R_ALIGN_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct RAlignStats : public AnalyzerStats
    {
        Sensitivity se;
    };

    struct RAlign
    {
        enum Level
        {
            Base,
            Exon,
            Splice,
        };

        struct Options : public AnalyzerOptions<RAlign::Level>
        {
            // Empty Implementation
        };

        static RAlignStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif