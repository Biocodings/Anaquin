#ifndef GI_D_ALIGN_HPP
#define GI_D_ALIGN_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct DAlignStats : public AnalyzerStats
    {
        Performance p;
        Counter c = DAnalyzer::counterSeqs();
    };

    struct DAlign
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static DAlignStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif