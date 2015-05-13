#ifndef GI_M_ASSEMBLY_HPP
#define GI_M_ASSEMBLY_HPP

#include "types.hpp"
#include "analyzer.hpp"

namespace Spike
{
    struct Contig
    {
        BasePair l;
        std::string id;
        std::string seq;
    };

    struct MAssemblyStats
    {
        std::vector<Contig> contigs;
        
        BasePair min, max;
        BasePair mean, sum;
        BasePair N20, N50, N80;
    };

    struct MAssembly
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static MAssemblyStats stats(const std::string &file);
        static MAssemblyStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif