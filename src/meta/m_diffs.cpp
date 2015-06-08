#include "meta/m_blast.hpp"
#include "meta/m_diffs.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MDiffs::Stats MDiffs::analyze(const std::string &file_1, const std::string &file_2, const Options &options)
{
    /*
     * The implementation is very similar to with one single sample. The only difference is that
     * we're interested in the log-fold ratio of the samples.
     */
    
    const auto stats_1 = Velvet::parse<MAssembly::Stats, Contig>(file_1);
    const auto stats_2 = Velvet::parse<MAssembly::Stats, Contig>(file_2);

    if (!options.psl_1.empty() && !options.psl_2.empty())
    {
        std::cout << "Using an aligment file: "  << options.psl_1 << std::endl;
        std::cout << "Using an aligment file: "  << options.psl_2 << std::endl;

        const auto r1 = MBlast::analyze(options.psl_1);
        const auto r2 = MBlast::analyze(options.psl_2);

        std::cout << "Creating an abundance plot" << std::endl;
        
        /*
         * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
         */
        
        std::vector<double> x, y;
        std::vector<std::string> z;
        
        for (const auto &meta : r1.metas)
        {
            // If the metaquin has an alignment
            if (!meta.second.aligns.empty())
            {
                // Known concentration
                const auto known = meta.second.seqB.abund() / meta.second.seqA.abund();
                
                BasePair measured = 0;

                /*
                 * Calculate measured concentration for this metaquin. The problem is that our
                 * alignment information is independent to the coverage. We'll need to link the
                 * pieces together. We'll also need to average out the contigs for the sequin.
                 */
                
                for (std::size_t i = 0; i < meta.second.temp.size(); i++)
                {
                    const auto &contig_1 = stats_1.contigs.at(meta.second.temp[i]);
                    const auto &contig_2 = stats_2.contigs.at(meta.second.temp[i]);

                    measured += (contig_2.k_cov / contig_2.seq.size()) / (contig_1.k_cov / contig_1.seq.size());
                    //measured += (contig_2.k_cov / meta.seqA.l.length()) / (contig_1.k_cov / meta.seqA.l.length());
                }

                x.push_back(log(known));
                y.push_back(log(measured));
                z.push_back(meta.second.id);
            }
        }
        
        // Generate a R script for a plot of abundance
        AnalyzeReporter::script("meta_diffs.R", x, y, z, options.writer);
    }

    return MDiffs::Stats();
}