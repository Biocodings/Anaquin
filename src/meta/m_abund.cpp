#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats MAbundance::analyze(const FileName &file, const MAbundance::Options &o)
{
    MAbundance::Stats stats;
    
    assert(!o.psl.empty());
    
    /*
     * Generate statistics for the alignment
     */
    
    o.info("Analyzing: " + o.psl);

    // Generate statistics for BLAT
    const auto bStats = MBlat::analyze(o.psl);

    o.info("Analyzing: " + file);

    // Generate statistics for Velvet, filtered by alignments (no use to keep non-synthetic in memory)
    const auto dnovo = Velvet::analyze<DAsssembly::Stats<Contig>, Contig>(file, &bStats);
    
    if (!dnovo.n)
    {
        throw std::runtime_error("No contig detected in the input file. Please check and try again.");
    }
    else if (dnovo.contigs.empty())
    {
        throw std::runtime_error("No contig aligned in the input file. Please check and try again.");
    }

    stats.n_chrT = dnovo.contigs.size();
    stats.n_hg38 = dnovo.n - stats.n_chrT;

    assert(stats.n_hg38 && stats.n_chrT);
    
    o.info("Analyzing the alignments");

    const auto &r = Standard::instance().r_meta;

    for (auto &meta : bStats.metas)
    {
        auto &align = meta.second;
        
        if (!r.match(align->seq->id))
        {
            o.warn((boost::format("%1% not defined in the mixture. Skipped.") % align->seq->id).str());
            continue;
        }
        
        /*
         * Calculate the limit of sensitivity. LOS is defined as the sequin with the lowest amount of
         * concentration while still detectable in the experiment.
         */

        if (stats.s.id.empty() || align->seq->abund(Mix_1, false) < stats.s.abund)
        {
            stats.s.id     = align->seq->id;
            stats.s.abund  = align->seq->abund(Mix_1, false);
            stats.s.counts = align->contigs.size();
        }
        
        /*
         * Plot the coverage relative to the known concentration for each assembled contig
         */
        
        if (!align->contigs.empty())
        {
            stats.h.at(meta.first)++;
            
            // Known concentration
            const auto known = align->seq->abund(Mix_1, false);
            
            /*
             * Measure concentration for this metaquin. Average out the coverage for each aligned contig.
             */

            Concentration measured = 0;
            
            for (auto i = 0; i < align->contigs.size(); i++)
            {
                if (!bStats.aligns.count(align->contigs[i].id))
                {
                    continue;
                }
                else if (!dnovo.contigs.count(align->contigs[i].id))
                {
                    o.warn((boost::format("%1% not found in the input file") % align->contigs[i].id).str());
                    continue;
                }
                
                const auto &contig = dnovo.contigs.at(align->contigs[i].id);

                assert(align->seq->l.length());
                assert(contig.k_cov && contig.k_len);
                
                switch (o.coverage)
                {
                    case KMerCov_Contig: { measured += contig.k_cov / contig.k_len;           break; }
                    case KMerCov_Sequin: { measured += contig.k_cov / align->seq->l.length(); break; }
                }

                /*
                 * Calculate for the average depth for alignment and sequin
                 */
                
                align->depthAlign  += align->contigs[i].l.length() * contig.k_cov / align->contigs[i].l.length();
                align->depthSequin += align->contigs[i].l.length() * contig.k_cov;
            }
            
            if (measured)
            {
                align->depthSequin = meta.second->depthSequin / align->seq->length;
                stats.add(align->seq->id, known, measured);
            }
        }
    }
    
    stats.ss = Standard::instance().r_meta.limit(stats.h);

    return stats;
}

void MAbundance::report(const FileName &file, const MAbundance::Options &o)
{
    const auto stats = MAbundance::analyze(file, o);

    o.info("Generating summary statistics");
    AnalyzeReporter::linear("MetaAbundance_summary.stats", stats, "contigs", o.writer, "sequins");
 
    o.info("Generating R script");
    AnalyzeReporter::scatter(stats,
                             "Expected abundance vs Measured coverage",
                             "MetaAbundance",
                             "Expected abudnance (attomol/ul)",
                             "Measured coverage (k-mer)",
                             "Expected abdunance (log2 attomol/ul)",
                             "Measured coverage (log2 k-mer)",
                             o.writer);
}