#include "meta/m_blast.hpp"
#include "meta/m_assembly.hpp"

using namespace Spike;

MAssembly::Stats MAssembly::analyze(const std::string &file, const Options &options)
{
    /*
     * The code for a specific assembler is indepenent to alignments. While it's
     * certinaly a good design, we'll need to link and bridge the details.
     */

    MAssembly::Stats stats;

    /*
     * Generate statistics from a specific assembler, references not required (de-novo assembly)
     */

    switch (options.tool)
    {
        case Velvet: { stats = Velvet::parse<MAssembly::Stats, Contig>(file); break; }
    }

    ModelStats ms;

    MBlast::Stats r; // TODO: Should this be here?

    if (!options.psl.empty())
    {
        throw std::invalid_argument("Alignment file needs to be specified");
    }

    std::cout << "Aligment file: " << options.psl << std::endl;
    
    // Analyse the given blast alignment file
    r = MBlast::analyze(options.psl);
    
    for (auto &meta : r.metas)
    {
        const auto &align = meta.second;
        
        // Ignore if there's a filter and the sequin is not one of those
        if (!options.filters.empty() && !options.filters.count(align.id))
        {
            continue;
        }
        
        /*
         * Calculate the limit of sensitivity. LOS is defined as the metaquin with the lowest amount of
         * concentration while still detectable in the experiment.
         */
        
        if (ms.s.id.empty() || align.seqA.abund() < ms.s.abund)
        {
            ms.s.id     = align.id;
            ms.s.abund  = align.seqA.abund();
            ms.s.counts = align.contigs.size();
        }
        
        /*
         * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
         */
        
        if (!align.contigs.empty())
        {
            // Known concentration
            const auto known = align.seqA.abund();
            
            /*
             * Calculate measured concentration for this metaquin. Average out the coverage for each aligned contig.
             */
            
            Concentration measured = 0;
            
            for (std::size_t i = 0; i < align.contigs.size(); i++)
            {
                // Crash if the alignment file doesn't match with the contigs...
                const auto &contig = stats.contigs.at(align.contigs[i].id);
                
                // Average relative to the size of the contig
                measured += contig.k_cov / contig.seq.size();
                
                // Average relative to the size of the sequin
                //measured += (double) contig.k_cov / meta.seqA.l.length();
                
                /*
                 * Calculate for the average depth for alignment and sequin
                 */
                
                meta.second.depthAlign  += align.contigs[i].l.length() * contig.k_cov / align.contigs[i].l.length();
                meta.second.depthSequin += align.contigs[i].l.length() * contig.k_cov;
            }
            
            meta.second.depthSequin = meta.second.depthSequin / align.seqA.length;
            assert(measured != 0);
            
            ms.z.push_back(align.id);
            ms.x.push_back(log(known));
            ms.y.push_back(log(measured));
        }
        
        assert(!ms.s.id.empty());
        
        // Generate a R script for a plot of abundance
        AnalyzeReporter::script("meta_abundance.R", ms.x, ms.y, ms.z, options.writer);
        
        std::cout << "Abundance plot generated" << std::endl;
    }

    /*
     * Write out results for each sequin
     */
    
    options.writer->open("abundance_stats.stats");
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

    options.writer->write((boost::format(format) % "ID"
                                                 % "Con"
                                                 % "Status"
                                                 % "DAlign"
                                                 % "DSequin"
                                                 % "Covered").str());
    
    for (auto &meta : r.metas)
    {
        const std::string status = meta.second.contigs.size() == 0 ? "Undetected" :
                                   meta.second.covered == 1.0      ? "Full" : "Partial";

        options.writer->write((boost::format(format) % meta.second.id
                                                     % meta.second.seqA.abund()
                                                     % status
                                                     % meta.second.depthAlign
                                                     % meta.second.depthSequin
                                                     % meta.second.covered).str());
    }
    
    options.writer->close();
    
    /*
     * Write out assembly results
     */

    options.writer->open("meta_assembly.stats");
    
    if (ms.x.size() <= 1)
    {
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
        
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
    }
    else
    {
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%";

        // Fit a simple linear regression model by maximum-likehihood
        const auto lm = ms.linear();

        options.writer->write((boost::format(format) % "Nodes"
                                                     % "N20"
                                                     % "N50"
                                                     % "N80"
                                                     % "min"
                                                     % "mean"
                                                     % "max"
                                                     % "total"
                                                     % "r"
                                                     % "slope"
                                                     % "r2"
                                                     % "los").str());
        options.writer->write((boost::format(format) % stats.contigs.size()
                                                     % stats.N20
                                                     % stats.N50
                                                     % stats.N80
                                                     % stats.min
                                                     % stats.mean
                                                     % stats.max
                                                     % stats.total
                                                     % lm.r
                                                     % lm.m
                                                     % lm.r2
                                                     % ms.s.abund).str());
    }
    
    options.writer->close();

    return stats;
}