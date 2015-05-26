#ifndef GI_M_ASSEMBLY_HPP
#define GI_M_ASSEMBLY_HPP

#include "analyzer.hpp"
#include "stats/denovo.hpp"

namespace Spike
{
    struct MAssemblyStats
    {
        DNStats dstats;
    };

    struct MAssembly
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static MAssemblyStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif