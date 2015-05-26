#include "stats/denovo.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MAssemblyStats MAssembly::analyze(const std::string &file, const Options &options)
{
    MAssemblyStats stats;

    stats.dstats = DNAsssembly::stats(file);
    
    
    return stats;;
}