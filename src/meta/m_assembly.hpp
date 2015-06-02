#ifndef GI_M_ASSEMBLY_HPP
#define GI_M_ASSEMBLY_HPP

#include "analyzer.hpp"
#include "stats/denovo.hpp"

namespace Spike
{
    struct MAssemblyStats
    {
        DNStats ds;
    };

    struct Node
    {
        std::string id;

        // The sequin that the node has been aligned to, empty if not aligned
        std::string sequin;

        // Coverage of the node
        double cov;
    };

    struct Velvet
    {
        struct VelvetStats
        {
            std::vector<Node> nodes;
        };

        static VelvetStats analyze(const std::string &contig, const std::string &blat);
    };
    
    struct MAssembly
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static MAssemblyStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif