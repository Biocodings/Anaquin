#include "meta/m_blast.hpp"
#include "meta/m_assembly.hpp"

using namespace Anaquin;

MAssembly::Stats MAssembly::analyze(const std::string &file, const Options &options)
{
    MAssembly::Stats stats;

    assert(!options.psl.empty());

    /*
     * Generate statistics for a specific assembler
     */

    switch (options.tool)
    {
        case Velvet: { stats = Velvet::parse<MAssembly::Stats, Contig>(file); break; }
    }

    LinearStats ms;

    // Analyse the given blast alignment file
    auto blat = MBlast::stats(options.psl);

    for (auto &meta : blat.metas)
    {
        const auto &align = meta.second;

        // Ignore if there's a filter and the sequin is not one of those
        if (!options.filters.empty() && !options.filters.count(align->seq->id))
        {
            continue;
        }
        
        /*
         * Calculate the limit of sensitivity. LOS is defined as the metaquin with the lowest amount of
         * concentration while still detectable in the experiment.
         */
        
        if (ms.s.id.empty() || align->seq->mixes.at(Mix_1) < ms.s.abund)
        {
            ms.s.id     = align->seq->id;
            ms.s.abund  = align->seq->mixes.at(Mix_1);
            ms.s.counts = align->contigs.size();
        }
        
        /*
         * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
         */
        
        if (!align->contigs.empty())
        {
            // Known concentration
            const auto known = align->seq->mixes.at(Mix_1);
            
            /*
             * Measure concentration for this metaquin. Average out the coverage for each aligned contig.
             */
            
            Concentration measured = 0;

            for (auto i = 0; i < align->contigs.size(); i++)
            {
                if (!stats.contigs.count(align->contigs[i].id))
                {
                    continue;
                }
                
                // Crash if the alignment file doesn't match with the contigs...
                const auto &contig = stats.contigs.at(align->contigs[i].id);

                // Average relative to the size of contig
                measured += contig.k_cov / contig.seq.size();
                
                // Average relative to the size of the sequin
                //measured += (double) contig.k_cov / meta.seqA.l.length();
                
                /*
                 * Calculate for the average depth for alignment and sequin
                 */
                
                meta.second->depthAlign  += align->contigs[i].l.length() * contig.k_cov / align->contigs[i].l.length();
                meta.second->depthSequin += align->contigs[i].l.length() * contig.k_cov;
            }
            
            if (measured)
            {
                meta.second->depthSequin = meta.second->depthSequin / align->seq->length;
                ms.add(align->seq->id, log2(known), log2(measured));
            }
        }
        
        assert(!ms.s.id.empty());
    }

    options.info("Generating linaer model");
    AnalyzeReporter::linear(ms, "MetaAssembly", "k-mer average", options.writer);
    
    /*
     * Write out results for each sequin
     */

    options.writer->open("MetaAssembly_quins.stats");
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

    options.writer->write((boost::format(format) % "ID"
                                                 % "con"
                                                 % "status"
                                                 % "avg_align"
                                                 % "avg_sequin"
                                                 % "covered").str());
    
    for (auto &meta : blat.metas)
    {
        const std::string status = meta.second->contigs.size() == 0 ? "Undetected" :
                                   meta.second->covered == 1.0      ? "Full" : "Partial";

        options.writer->write((boost::format(format) % meta.second->seq->id
                                                     % meta.second->seq->mixes.at(Mix_1)
                                                     % status
                                                     % meta.second->depthAlign
                                                     % meta.second->depthSequin
                                                     % meta.second->covered).str());
    }
    
    options.writer->close();

    {
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

        options.writer->open("MetaAssembly_summary.stats");
        options.writer->write((boost::format(format) % "contigs"
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
    }

    return stats;
}