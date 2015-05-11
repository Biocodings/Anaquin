#ifndef GI_D_ALIGN_HPP
#define GI_D_ALIGN_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct DAlignStats
    {
        // Overall performance for the exon-level
        Confusion me;

        // Overlal LOS for the exon-level
        Sensitivity se;

        Counter ce = DAnalyzer::sequinCounter();
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