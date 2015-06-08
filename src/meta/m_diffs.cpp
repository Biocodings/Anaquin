#include "meta/m_blast.hpp"
#include "meta/m_diffs.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MDiffs::Stats MDiffs::analyze(const std::string &file_1, const std::string &file_2, const Options &options)
{
    /*
     * The implementation is very similar to one single sample. The only difference is that
     * we're interested in the log-fold change of the samples.
     */
    
    const auto stats_1 = Velvet::parse<MAssembly::Stats, Contig>(file_1);
    const auto stats_2 = Velvet::parse<MAssembly::Stats, Contig>(file_2);

    if (!options.psl_1.empty() && !options.psl_2.empty())
    {
        std::cout << "Using an aligment file: "  << options.psl_1 << std::endl;
        std::cout << "Using an aligment file: "  << options.psl_2 << std::endl;

        const auto r1 = MBlast::analyze(options.psl_1);
        const auto r2 = MBlast::analyze(options.psl_2);

        std::cout << "Creating a differential plot" << std::endl;
        
        /*
         * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
         */
        
        std::vector<double> x, y;
        std::vector<std::string> z;
        
        // Marginal observation for mixture A
        std::map<SequinID, Concentration> x1;
        
        // Marginal observation for mixture B
        std::map<SequinID, Concentration> x2;

        for (const auto &meta : r1.metas)
        {
            const auto &align = meta.second;
            
            // If the metaquin has an alignment
            if (!align.contigs.empty())
            {
                /*
                 * Calculate measured concentration for this metaquin. Average out
                 * the coverage for each aligned contig.
                 */
                
                Concentration measured = 0;
                
                for (std::size_t i = 0; i < align.contigs.size(); i++)
                {
                    const auto &contig = stats_1.contigs.at(align.contigs[i].id);

                    // Average relative to the size of the contig
                    measured += contig.k_cov / contig.seq.size();
                    
                    // Average relative to the size of the sequin
                    //measured += (double) contig.k_cov / meta.seqA.l.length();
                }
                
                assert(measured != 0);
                x1[align.id] = measured;
            }
            else
            {
                x1[align.id] = 0;
            }
        }

        for (const auto &meta : r2.metas)
        {
            const auto &align = meta.second;
            
            // If the metaquin has an alignment
            if (!align.contigs.empty())
            {
                /*
                 * Calculate measured concentration for this metaquin. Average out
                 * the coverage for each aligned contig.
                 */
                
                Concentration measured = 0;
                
                for (std::size_t i = 0; i < align.contigs.size(); i++)
                {
                    const auto &contig = stats_2.contigs.at(align.contigs[i].id);
                    
                    // Average relative to the size of the contig
                    measured += contig.k_cov / contig.seq.size();
                    
                    // Average relative to the size of the sequin
                    //measured += (double) contig.k_cov / meta.seqA.l.length();
                }
                
                assert(measured != 0);
                x2[align.id] = measured;
            }
            else
            {
                x2[align.id] = 0;
            }
        }
        
        assert(x1.size() == x2.size());
        assert(r1.metas.size() == r2.metas.size());
        
        for (const auto &meta : r1.metas)
        {
            const auto &align = meta.second;

            // If the metaquin has an alignment
            if (!align.contigs.empty())
            {
                // Known concentration
                const auto known = align.seqB.abund() / align.seqA.abund();
                
                // TODO: ????
                if (x2.at(align.id) && x1.at(align.id))
                {
                    // Ratio of the marginal concentration
                    const auto measured = x2.at(align.id) / x1.at(align.id);
                    
                    x.push_back(log(known));
                    y.push_back(log(measured));
                    z.push_back(align.id);                    
                }
            }
        }
        
        // Generate a R script for a plot of abundance
        AnalyzeReporter::script("meta_diffs.R", x, y, z, options.writer);
    }

    return MDiffs::Stats();
}