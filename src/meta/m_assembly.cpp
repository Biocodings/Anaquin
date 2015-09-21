#include "meta/m_blat.hpp"
#include "meta/m_assembly.hpp"

using namespace Anaquin;

MAssembly::Stats MAssembly::analyze(const FileName &file, const Options &o)
{
    MAssembly::Stats stats;

    assert(!o.psl.empty());

    /*
     * Generate statistics for the assembler
     */

    o.info("Parsing the contig file: " + file);

    switch (o.tool)
    {
        case Velvet: { stats = Velvet::parse<MAssembly::Stats, Contig>(file); break; }
    }

    o.info("Analyzing the PSL file");
    
    // Analyse the blat alignment file
    stats.blat = MBlat::analyze(o.psl);

    o.info("Analyzing the PSL alignments");

    for (auto &meta : stats.blat.metas)
    {
        auto &align = meta.second;

        // Ignore if there's a filter and the sequin is not one of those
        if (!o.filters.empty() && !o.filters.count(align->seq->id))
        {
            continue;
        }
        
        /*
         * Calculate the limit of sensitivity. LOS is defined as the sequin with the lowest amount of
         * concentration while still detectable in the experiment.
         */
        
        const auto needMixture = align->seq->mixes.count(Mix_1);
        
        if (needMixture)
        {
            if (stats.lm.s.id.empty() || align->seq->mixes.at(Mix_1) < stats.lm.s.abund)
            {
                stats.lm.s.id     = align->seq->id;
                stats.lm.s.abund  = align->seq->mixes.at(Mix_1);
                stats.lm.s.counts = align->contigs.size();
            }
        }
        
        /*
         * Plot the coverage relative to the known concentration for each assembled contig
         */

        if (needMixture && !align->contigs.empty())
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
                measured += contig.k_cov / contig.k_len;
                
                // Average relative to the size of the sequin
                //measured += (double) contig.k_cov / meta.seqA.l.length();
                
                /*
                 * Calculate for the average depth for alignment and sequin
                 */

                align->depthAlign  += align->contigs[i].l.length() * contig.k_cov / align->contigs[i].l.length();
                align->depthSequin += align->contigs[i].l.length() * contig.k_cov;
            }
            
            if (measured)
            {
                align->depthSequin = meta.second->depthSequin / align->seq->length;
                stats.lm.add(align->seq->id, log2(known), log2(measured));
            }
        }
    }
    
    return stats;
}

MAssembly::Stats MAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = MAssembly::analyze(file, o);

    o.logInfo("Generating summary statistics");
    
    /*
     * Generate summary statistics
     */

    {
        o.writer->open("MetaAssembly_summary.stats");
        
        const auto summary = "Summary for dataset: %1%\n\n"
                             "   Community: %2%\n"
                             "   Synthetic: %3%\n\n"
                             "   Contigs: %4%\n"
                             "   Assembled: %5%\n"
                             "   Reference: %6%\n\n"
                             "   ***\n"
                             "   *** The following statistics on the synthetic community\n"
                             "   ***\n\n"
                             "   Contigs:  %7%\n"
                             "   N20:    %8%\n"
                             "   N50: %9%\n"
                             "   N80: %10%\n"
                             "   min: %11%\n"
                             "   mean: %12%\n"
                             "   max: %13%\n"
                             "   ***\n"
                             "   *** The following overlapping statistics are computed by proportion\n"
                             "   ***\n\n"
                             "   Match: %14%\n"
                             "   Gaps: %15%\n"
                             "   Mismatch: %16%\n";
        
        o.writer->write((boost::format(summary) % file
                                                % stats.blat.n_hg38
                                                % stats.blat.n_chrT
                                                % stats.blat.aligns.size()
                                                % stats.blat.sequin()
                                                % stats.blat.metas.size()
                                                % stats.contigs.size()
                                                % stats.N20
                                                % stats.N50
                                                % stats.N80
                                                % stats.min
                                                % stats.mean
                                                % stats.max
                                                % stats.blat.overMatch()
                                                % stats.blat.overGaps()
                                                % stats.blat.overMismatch()).str());
        o.writer->close();
    }
    
    /*
     * Generate results for each sequin
     */

    o.logInfo("Generating sequins statistics");

    {
        o.writer->open("MetaAssembly_quins.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.writer->write((boost::format(format) % "id"
                                               % "contigs"
                                               % "covered"
                                               % "mismatch"
                                               % "gaps").str());
        for (const auto &i : stats.blat.metas)
        {
            const auto &align = i.second;
            
            o.writer->write((boost::format(format) % align->seq->id
                                                   % align->contigs.size()
                                                   % align->covered
                                                   % align->mismatch
                                                   % align->gaps).str());
        }
        
        o.writer->close();
    }

    return stats;
}