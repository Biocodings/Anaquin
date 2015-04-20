#ifndef GI_DALIGNER_HPP
#define GI_DALIGNER_HPP

#include "analyzer.hpp"
#include "classify.hpp"

namespace Spike
{
    struct DAlignerStats : public AnalyzerStats
    {
        // Empty Implementation
    };

    struct DAligner
    {
        enum Mode
        {
            Base,
            Exon,
            Splice,
        };

        struct Options : public AnalyzerOptions<DAligner::Mode>
        {
            // Empty Implementation
        };

        static DAlignerStats analyze(const std::string &file, const DAligner::Options &options = DAligner::Options());
    };
}

#endif