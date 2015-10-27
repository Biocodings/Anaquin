#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats MAbundance::analyze(const FileName &file, const MAbundance::Options &o)
{
    MAbundance::logOptions(o);
    
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
    const auto dStats = Velvet::analyze<DAsssembly::Stats<Contig>, Contig>(file, &bStats);
    
    if (!dStats.n)
    {
        throw std::runtime_error("No contig detected in the input file. Please check and try again.");
    }
    else if (dStats.contigs.empty())
    {
        throw std::runtime_error("No contig aligned in the input file. Please check and try again.");
    }

    stats.n_chrT = dStats.contigs.size();
    stats.n_expT = dStats.n - stats.n_chrT;

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
        
        const auto p = MAbundance::calculate(stats, bStats, dStats, align->seq->id, *meta.second, o, o.coverage);
        
        if (p.x && p.y)
        {
            stats.add(align->seq->id, p.x, p.y);
        }
    }
    
    stats.ss = Standard::instance().r_meta.limit(stats.h);

    return stats;
}

void MAbundance::report(const FileName &file, const MAbundance::Options &o)
{
    const auto stats = MAbundance::analyze(file, o);

    o.info("Generating summary statistics");
    AnalyzeReporter::linear("MetaAbundance_summary.stats", file, stats, "contigs", o.writer, "sequins", "Community");
 
    o.info("Generating Bioconductor");
    AnalyzeReporter::scatter(stats,
                             "Expected abundance vs Measured coverage",
                             "MetaAbundance",
                             "Expected abudnance (attomol/ul)",
                             "Measured coverage (k-mer)",
                             "Expected abdunance (log2 attomol/ul)",
                             "Measured coverage (log2 k-mer)",
                             o.writer);
}