#include "meta/m_blast.hpp"
#include "meta/m_assembly.hpp"

using namespace Anaquin;

MAssembly::Stats MAssembly::analyze(const std::string &file, const Options &o)
{
    MAssembly::Stats stats;

    assert(!o.psl.empty());

    /*
     * Generate statistics for the assembler
     */

    switch (o.tool)
    {
        case Velvet: { stats = Velvet::parse<MAssembly::Stats, Contig>(file); break; }
    }

    // Analyse the blat alignment file
    stats.blat = MBlast::stats(o.psl);

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

        if (stats.lm.s.id.empty() || align->seq->mixes.at(Mix_1) < stats.lm.s.abund)
        {
            stats.lm.s.id     = align->seq->id;
            stats.lm.s.abund  = align->seq->mixes.at(Mix_1);
            stats.lm.s.counts = align->contigs.size();
        }
        
        /*
         * Plot the coverage relative to the known concentration for each assembled contig
         */

        if (!align->contigs.empty())
        {
            // Known concentration (TODO: We shouldn't use the data structure...)
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

        assert(!stats.lm.s.id.empty());
    }
    
    return stats;
}

MAssembly::Stats MAssembly::report(const std::string &file, const Options &o)
{
    const auto stats = MAssembly::analyze(file, o);

    /*
     * Generate summary statistics
     */

    {
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

        o.writer->open("MetaAssembly_summary.stats");
        
        const auto summary = "Summary for dataset: %1%\n\n"
                             "   Number of contigs:  %2%\n"
                             "   N20:    %3%\n"
                             "   N50: %4%\n"
                             "   N80: %5%\n"
                             "   min: %6%\n"
                             "   mean: %7%\n"
                             "   max: %8%\n"
                             "   total: %9%\n"
        ;
        
        o.writer->write((boost::format(summary) % file
                                                % stats.contigs.size()
                                                % stats.N20
                                                % stats.N50
                                                % stats.N80
                                                % stats.min
                                                % stats.mean
                                                % stats.max
                                                % stats.total).str());
        o.writer->close();
    }

    /*
     * Generate results for each sequin
     */

    {
        o.writer->open("MetaAssembly_quins.stats");
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

        o.writer->write((boost::format(format) % "ID"
                                               % "con"
                                               % "status"
                                               % "avg_align"
                                               % "avg_sequin"
                                               % "covered").str());
        
        for (const auto &meta : stats.blat.metas)
        {
            const auto &align = meta.second;
            const auto detect = align->contigs.size() != 0;
            const auto status = detect ? std::to_string(align->covered) : "-";

            o.writer->write((boost::format(format) % align->seq->id
                                                   % align->seq->mixes.at(Mix_1)
                                                   % status
                                                   % (detect ? std::to_string(align->depthAlign)  : "-")
                                                   % (detect ? std::to_string(align->depthSequin) : "-")
                                                   % (detect ? std::to_string(align->covered)     : "-")
                             ).str());
        }
    }

    o.writer->close();

    return stats;
}