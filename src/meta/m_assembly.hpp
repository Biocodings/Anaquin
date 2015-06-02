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

    struct Velvet
    {
        struct VelvetStats
        {
            
        };

        static VelvetStats analyze(const std::string &file);
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