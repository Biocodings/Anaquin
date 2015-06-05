#include "meta/m_blast.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MAssembly::Stats MAssembly::analyze(const std::string &file, const Options &options)
{
    MAssembly::Stats stats = Velvet::parse<MAssembly::Stats, MAssembly::AlignedContig>(file);

    // Prefer the alignment file from user if provided
    if (!options.blast.empty())
    {
        std::cout << "Using an aligment file: "  << options.blast << std::endl;
        
        // Analyse the given blast alignment file
        const auto r = MBlast::analyze(options.blast);
        
        std::cout << "Creating a abundance plot" << std::endl;

        /*
         * Plot the coverage realtive to the known concentration (in attamoles/ul) of each assembled contig.
         */

        std::vector<double> x, y;
        std::vector<std::string> z;

        for (const auto &contig : stats.contigs)
        {
            // If the contig has an alignment
            if (r.aligns.count(contig.id))
            {
                const auto &seq = r.aligns.at(contig.id).seqA;
                
                // Known concentration from the matching sequin
                const auto known = seq.abund();

                assert(contig.k_cov);
                
                // Measured coverage (alignment coverage)
                const auto measured = contig.k_cov / seq.l.length();
                //const auto measured = contig.k_cov contig.seq.length();
                
                x.push_back(log(known));
                y.push_back(log(measured));
                z.push_back(contig.id);
            }
        }

        // Generate a R script for a plot of abundance
        AnalyzeReporter::script("abundance.R", x, y, z, options.writer);
    }
    
    /*
     * Write out assembly results
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";

    options.writer->open("assembly.stats");
    options.writer->write((boost::format(format) % "Nodes" % "N20" % "N50" % "N80" % "min" % "mean" % "max" % "total" % "reads").str());
    options.writer->write((boost::format(format) % stats.contigs.size()
                                                 % stats.N20
                                                 % stats.N50
                                                 % stats.N80
                                                 % stats.min
                                                 % stats.mean
                                                 % stats.max
                                                 % stats.total
                                                 % -1).str());
    options.writer->close();

    return stats;
}