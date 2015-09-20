#include "meta/m_diffs.hpp"
#include "meta/m_assembly.hpp"

using namespace Anaquin;

MDiffs::Stats MDiffs::report(const FileName &file_1, const FileName &file_2, const Options &o)
{
    MDiffs::Stats stats;

    /*
     * The implementation is very similar to one single sample. The only difference is that
     * we're interested in the log-fold change of the samples.
     */
    
    const auto stats_1 = Velvet::parse<MAssembly::Stats, Contig>(file_1);
    const auto stats_2 = Velvet::parse<MAssembly::Stats, Contig>(file_2);

    if (!o.pA.empty() && !o.pB.empty())
    {
        o.info((boost::format("Using alignment: %1%") % o.pA).str());
        o.info((boost::format("Using alignment: %1%") % o.pB).str());

        const auto r1 = MBlast::analyze(o.pA);
        const auto r2 = MBlast::analyze(o.pB);

        o.info("Creating a differential plot");
        
        /*
         * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
         */
        
        // Marginal for mixture A
        std::map<SequinID, Coverage> y1;
        
        // Marginal for mixture B
        std::map<SequinID, Coverage> y2;

        for (const auto &meta : r1.metas)
        {
            const auto &align = meta.second;
            
            // If the metaquin has an alignment
            if (!align->contigs.empty())
            {
                /*
                 * Calculate measured concentration for this metaquin. Average out
                 * the coverage for each aligned contig.
                 */
                
                Concentration measured = 0;
                
                for (std::size_t i = 0; i < align->contigs.size(); i++)
                {
                    const auto &contig = stats_1.contigs.at(align->contigs[i].id);

                    // Average relative to the size of the contig
                    measured += contig.k_cov / contig.k_len;
                    
                    // Average relative to the size of the sequin
                    //measured += contig.k_cov / meta.second.seqA.length;
                }
                
                assert(measured != 0);
                y1[align->id()] = measured;
            }
            else
            {
                y1[align->id()] = 0;
            }
        }

        for (const auto &meta : r2.metas)
        {
            const auto &align = meta.second;
            
            // If the metaquin has an alignment
            if (!align->contigs.empty())
            {
                /*
                 * Calculate measured concentration for this metaquin. Average out
                 * the coverage for each aligned contig.
                 */
                
                Concentration measured = 0;
                
                for (std::size_t i = 0; i < align->contigs.size(); i++)
                {
                    const auto &contig = stats_2.contigs.at(align->contigs[i].id);
                    
                    // Average relative to the size of the contig
                    measured += contig.k_cov / contig.k_len;
                    
                    // Average relative to the size of the sequin
                    //measured += contig.k_cov / meta.second.seqB.length;
                }
                
                assert(measured != 0);
                y2[align->id()] = measured;
            }
            else
            {
                y2[align->id()] = 0;
            }
        }
        
        assert(y1.size() == y2.size());
        assert(r1.metas.size() == r2.metas.size());
        
        for (const auto &meta : r1.metas)
        {
            const auto &align = meta.second;

            // If the metaquin has an alignment
            if (!align->contigs.empty())
            {
                // Only when the sequin has mapping for both mixtures...
                if (y2.at(align->id()) && y1.at(align->id()))
                {
                    // Ignore if there's a filter and the sequin is not one of those
                    if (!o.filters.empty() && !o.filters.count(align->id()))
                    {
                        continue;
                    }
                    
                    // Known concentration
                    const auto known = align->seq->mixes.at(Mix_2) / align->seq->mixes.at(Mix_1);

                    // Ratio of the marginal concentration
                    const auto measured = y2.at(align->id()) / y1.at(align->id());
                    
                    stats.add(align->id(), log2(known), log2(measured));

                    SequinDiff d;
                    
                    d.id   = align->id();
                    d.ex_A = align->seq->mixes.at(Mix_1);
                    d.ex_B = align->seq->mixes.at(Mix_2);
                    d.ob_A = y1.at(align->id());
                    d.ob_B = y2.at(align->id());
                    d.ex_fold = known;
                    d.ob_fold = measured;

                    stats.diffs.push_back(d);
                }
            }
        }

        //AnalyzeReporter::linear(stats, "m_diffs", "k-mer average", options.writer);
    }
    
    o.info("Generating statistics");
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";

    o.writer->open("diff.stats");
    o.writer->write((boost::format(format) % "id"
                                           % "expect_A"
                                           % "expect_B"
                                           % "measure_A"
                                           % "measure_B"
                                           % "expect_fold"
                                           % "measure_fold").str());

    for (const auto &diff : stats.diffs)
    {
        o.writer->write((boost::format(format) % diff.id
                                               % diff.ex_A
                                               % diff.ob_A
                                               % diff.ex_B
                                               % diff.ob_B
                                               % diff.ex_fold
                                               % diff.ob_fold).str());
    }
    
    o.writer->close();
    
    return stats;
}