#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats MAbundance::analyze(const FileName &file, const MAbundance::Options &o)
{
    MAbundance::Stats stats;
    
    assert(!o.psl.empty());
    
    /*
     * Generate statistics for the alignment
     */
    
    o.info("Analyzing PSL: " + o.psl);

    // Generate statistics for BLAT
    const auto bStats = MBlat::analyze(o.psl);

    o.info("Analyzing contig: " + file);

    // Generate statistics for Velvet
    const auto vStats = Velvet::analyze<DAsssembly::Stats<Contig>, Contig>(file);
    
    o.info("Analyzing the PSL alignments");
    
    for (auto &meta : bStats.metas)
    {
        auto &align = meta.second;
        
        /*
         * Calculate the limit of sensitivity. LOS is defined as the sequin with the lowest amount of
         * concentration while still detectable in the experiment.
         */
        
        if (stats.s.id.empty() || align->seq->mixes.at(Mix_1) < stats.s.abund)
        {
            stats.s.id     = align->seq->id;
            stats.s.abund  = align->seq->mixes.at(Mix_1);
            stats.s.counts = align->contigs.size();
        }
        
        /*
         * Plot the coverage relative to the known concentration for each assembled contig
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
                if (!bStats.aligns.count(align->contigs[i].id))
                {
                    continue;
                }
                
                // Crash if the alignment file doesn't match with the contigs...
                const auto &contig = vStats.contigs.at(align->contigs[i].id);
                
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
    
    return stats;
}

void MAbundance::report(const FileName &file, const MAbundance::Options &o)
{
    const auto stats = MAbundance::analyze(file, o);

    o.info("Generating linaer model");
    AnalyzeReporter::linear("MetaAbundance_summary.stats", stats, "...", o.writer);
    
    o.info("Generating linaer model");
    AnalyzeReporter::scatter(stats, "MetaAbundance", "Expected abudnance (attomol/ul)", "K-Mer average", o.writer);
}