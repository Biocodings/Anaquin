#include "tokens.hpp"
#include "stats/denovo.hpp"
#include "meta/m_assembly.hpp"
#include "parsers/parser_fa.hpp"

using namespace Spike;

Velvet::VelvetStats Velvet::analyze(const std::string &file)
{
    /*
     * Read coverage from the contig file. The format would be like:
     *
     *      >NODE_77460_length_31_cov_1.129032
     */

    std::vector<std::string> tokens;
    
    ParserFA::parse(file, [&](const FALine &l, const ParserProgress &)
    {
        Tokens::split(l.id, "_", tokens);

        const auto cov = stof(tokens[tokens.size() - 1]);
        
        
        
        std::cout << l.id << std::endl;
    });

    
    
    
    
    return Velvet::VelvetStats();
}

MAssemblyStats MAssembly::analyze(const std::string &file, const Options &options)
{
    MAssemblyStats stats;

    // Calculate the general statistics for de-novo assembly
    stats.ds = DNAsssembly::stats(file);

    /*
     * Write out assembly results
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

    options.writer->open("meta_assembly.stats");
    options.writer->write((boost::format(format) % "N20" % "N50" % "N80" % "min" % "mean" % "max").str());
    options.writer->write((boost::format(format) % stats.ds.N20
                                                 % stats.ds.N50
                                                 % stats.ds.N80
                                                 % stats.ds.min
                                                 % stats.ds.mean
                                                 % stats.ds.max).str());
    options.writer->close();

    return stats;;
}