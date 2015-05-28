#include "stats/denovo.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MAssemblyStats MAssembly::analyze(const std::string &file, const Options &options)
{
    MAssemblyStats stats;

    // Calculate the general statistics for de-novo assembly
    stats.dstats = DNAsssembly::stats(file);
    
    /*
     * Write out assembly results
     */
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

    options.writer->open("meta_assembly.stats");
    options.writer->write((boost::format(format) % "N20" % "N50" % "N80" % "min" % "mean" % "max").str());
    options.writer->write((boost::format(format) % stats.dstats.N20
                                                 % stats.dstats.N50
                                                 % stats.dstats.N80
                                                 % stats.dstats.min
                                                 % stats.dstats.mean
                                                 % stats.dstats.max).str());
    options.writer->close();

    return stats;;
}