#include "meta/m_diffs.hpp"
#include "meta/m_assembly.hpp"

using namespace Anaquin;

MDiffs::Stats MDiffs::report(const FileName &file_1, const FileName &file_2, const Options &o)
{
    MDiffs::Stats stats;

    assert(!o.pA.empty() && !o.pB.empty());
    
    o.info((boost::format("Analyzing: %1%") % o.pA).str());
    stats.align_1 = MBlat::analyze(o.pA);

    o.info((boost::format("Analyzing: %1%") % o.pB).str());
    stats.align_2 = MBlat::analyze(o.pB);

    /*
     * The implementation is very similar to a single sample. The only difference is that
     * we're interested in the log-fold change.
     */
    
    o.info((boost::format("Analyzing: %1%") % file_1).str());
    const auto stats_1 = Velvet::analyze<MAssembly::Stats, Contig>(file_1, &stats.align_1);

    o.info((boost::format("Analyzing: %1%") % file_2).str());
    const auto stats_2 = Velvet::analyze<MAssembly::Stats, Contig>(file_2, &stats.align_2);

    /*
     * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
     */
    
    // Marginal for mixture A
    std::map<SequinID, Coverage> y1;
    
    // Marginal for mixture B
    std::map<SequinID, Coverage> y2;
    
    /*
     * Analyzing for the first sample
     */
    
    o.info("Analyzing for the first sample");

    for (const auto &meta : stats.align_1.metas)
    {
        const auto &align = meta.second;
        
        // If there is an alignment for the sequin...
        if (!align->contigs.empty())
        {
            /*
             * Average out the coverage for each aligned contig.
             */
            
            Concentration measured = 0;
            
            for (auto i = 0; i < align->contigs.size(); i++)
            {
                const auto &contig = stats_1.contigs.at(align->contigs[i].id);

                switch (o.coverage)
                {
                    case MAbundance::KMerCov_Contig: { measured += contig.k_cov / contig.k_len;           break; }
                    case MAbundance::KMerCov_Sequin: { measured += contig.k_cov / align->seq->l.length(); break; }
                }
            }
            
            assert(measured != 0);
            y1[align->id()] = measured;
        }
        else
        {
            y1[align->id()] = 0;
        }
    }

    /*
     * Analyzing for the second sample
     */

    o.info("Analyzing for the second sample");
    
    for (const auto &meta : stats.align_2.metas)
    {
        const auto &align = meta.second;
        
        // If there is an alignment for the sequin...
        if (!align->contigs.empty())
        {
            Concentration measured = 0;
            
            for (auto i = 0; i < align->contigs.size(); i++)
            {
                const auto &contig = stats_2.contigs.at(align->contigs[i].id);
                
                switch (o.coverage)
                {
                    case MAbundance::KMerCov_Contig: { measured += contig.k_cov / contig.k_len;           break; }
                    case MAbundance::KMerCov_Sequin: { measured += contig.k_cov / align->seq->l.length(); break; }
                }
            }
            
            assert(measured != 0);
            y2[align->id()] = measured;
        }
        else
        {
            y2[align->id()] = 0;
        }
    }

    /*
     * Merging information, note that we can only do a differential comparison if the sequin appears
     * detected in both samples
     */
    
    for (const auto &meta : stats.align_1.metas)
    {
        const auto &align = meta.second;
        
        // If the metaquin has an alignment
        if (!align->contigs.empty())
        {
            // Only when the align is detected in both samples
            if (y2.at(align->id()) && y1.at(align->id()))
            {
                // Known concentration
                const auto known = align->seq->abund(Mix_2) / align->seq->abund(Mix_1);

                // Ratio of the marginal concentration
                const auto measured = y2.at(align->id()) / y1.at(align->id());
                
                stats.add(align->id(), known, measured);
                
                SequinDiff d;
                
                d.id   = align->id();
                d.ex_A = align->seq->abund(Mix_1);
                d.ex_B = align->seq->abund(Mix_2);
                d.ob_A = y1.at(align->id());
                d.ob_B = y2.at(align->id());
                d.ex_fold = known;
                d.ob_fold = measured;

                stats.diffs.insert(d);
            }
        }
    }

    o.info((boost::format("Detected %1% sequins in estimating differential") % stats.size()).str());

    /*
     * Generating summary statistics
     */
    
    o.info("Generating statistics");

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";

    o.writer->open("MetaDifferent_summary.stats");
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
    
    /*
     * Generating detailed statistics for each sequin
     */
    
    o.writer->open("MetaDifferent_quin.stats");
    
    o.writer->close();
    
    /*
     * Generating differential comparisons for both samples
     */

    o.info("Generating summary statistics");
    AnalyzeReporter::linear("MetaDifferent_summary.stats", stats, "sequins", o.writer);

    /*
     * Generate an R script
     */
    
    o.info("Generating R script");
    AnalyzeReporter::scatter(stats, "MetaDifferent", "Expected log2 fold change (log2)", "Measured log2 fold change", o.writer);
    
    return stats;
}