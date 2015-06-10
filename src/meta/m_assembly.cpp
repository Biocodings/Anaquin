#include "meta/m_blast.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MAssembly::Stats MAssembly::analyze(const std::string &file, const Options &options)
{
    /*
     * The code for a specific assembler is indepenent to alignment for contigs.
     * While it is certinaly a good design, we'll need to link the information.
     */

    MAssembly::Stats stats;

    /*
     * Generate statistics for contigs, no reference sequins needed
     */

    switch (options.tool)
    {
        case Velvet: { stats = Velvet::parse<MAssembly::Stats, Contig>(file); break; }
    }

    // Prefer a user supplied alignment file if any
    if (!options.psl.empty())
    {
        std::cout << "Using an aligment file: " << options.psl << std::endl;

        // Analyse the given blast alignment file
        const auto r = MBlast::analyze(options.psl);

        std::cout << "Creating an abundance plot" << std::endl;

        /*
         * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
         */

        std::vector<double> x, y;
        std::vector<std::string> z;

        for (const auto &meta : r.metas)
        {
            const auto &align = meta.second;

            // Ignore if there's a filter and the sequin is not one of those
            if (!options.filters.empty() && !options.filters.count(align.id))
            {
                continue;
            }

            // If the metaquin has an alignment
            if (!align.contigs.empty())
            {
                // Known concentration
                const auto known = align.seqA.abund();

                /*
                 * Calculate measured concentration for this metaquin. Average out
                 * the coverage for each aligned contig.
                 */

                Concentration measured = 0;

                for (std::size_t i = 0; i < align.contigs.size(); i++)
                {
                    // Crash if the alignment file doesn't match with the contigs
                    const auto &contig = stats.contigs.at(align.contigs[i].id);
                    
                    // Average relative to the size of the contig
                    measured += contig.k_cov / contig.seq.size();
                    
                    // Average relative to the size of the sequin
                    //measured += (double) contig.k_cov / meta.seqA.l.length();
                }
                
                assert(measured != 0);

                x.push_back(log(known));
                y.push_back(log(measured));
                z.push_back(align.id);
            }
        }

        // Generate a R script for a plot of abundance
        AnalyzeReporter::script("meta_abundance.R", x, y, z, options.writer);
        
        std::cout << "Abundance plot generated" << std::endl;
    }

    /*
     * Write out assembly results
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

    options.writer->open("assembly.stats");
    options.writer->write((boost::format(format) % "Nodes"
                                                 % "N20"
                                                 % "N50"
                                                 % "N80"
                                                 % "min"
                                                 % "mean"
                                                 % "max"
                                                 % "total").str());
    options.writer->write((boost::format(format) % stats.contigs.size()
                                                 % stats.N20
                                                 % stats.N50
                                                 % stats.N80
                                                 % stats.min
                                                 % stats.mean
                                                 % stats.max
                                                 % stats.total).str());
    options.writer->close();

    return stats;
}