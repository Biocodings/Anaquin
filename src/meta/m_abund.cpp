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
    auto t = MBlat::analyze(o.psl);

    /*
     * Generate statistics for the assembly
     */
 
    o.info("Analyzing: " + file);

    switch (o.tool)
    {
        case Velvet:  { stats.assembly = Velvet::analyze<MAssembly::Stats, Contig>(file, &t);             break; }
        case RayMeta: { stats.assembly = RayMeta::analyze<MAssembly::Stats, Contig>(file, o.contigs, &t); break; }
    }

    stats.blat = t;

    if (!stats.assembly.n)
    {
        throw std::runtime_error("No contig detected in the input file. Please check and try again.");
    }
    else if (stats.assembly.contigs.empty())
    {
        throw std::runtime_error("No contig aligned in the input file. Please check and try again.");
    }

    stats.n_chrT = stats.assembly.contigs.size();
    stats.n_expT = stats.assembly.n - stats.n_chrT;

    o.info("Analyzing the alignments");

    const auto &r = Standard::instance().r_meta;

    for (auto &meta : stats.blat.metas)
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
        
        const auto p = MAbundance::calculate(stats, stats.blat, stats.assembly, align->seq->id, *meta.second, o, o.coverage);
        
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
    AnalyzeReporter::linear("MetaAbund_summary.stats", file, stats, "contigs", o.writer, "sequins");
 
    o.info("Generating Bioconductor");
    AnalyzeReporter::scatter(stats,
                             "Expected abundance vs Measured coverage",
                             "MetaAbundance",
                             "Expected abudnance (attomol/ul)",
                             "Measured coverage (k-mer)",
                             "Expected abdunance (log2 attomol/ul)",
                             "Measured coverage (log2 k-mer)",
                             o.writer);
    
    /*
     * Generating detailed statistics for each contig
     */
    
    {
        o.writer->open("MetaAbund_contigs.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";

        o.writer->write((boost::format(format) % "ID"
                                               % "Sequin ID"
                                               % "Length"
                                               % "Coverage"
                                               % "Normalized").str());
        
        for (const auto &i : stats.blat.aligns)
        {
            if (stats.assembly.contigs.count(i.first))
            {
                const auto &contig = stats.assembly.contigs.at(i.first);
                
                o.writer->write((boost::format(format) % i.first
                                                       % i.second->id()
                                                       % contig.k_len
                                                       % contig.k_cov
                                                       % contig.normalized()
                                 ).str());
            }
            else
            {
                o.writer->write((boost::format(format) % i.first
                                                       % i.second->id()
                                                       % "-"
                                                       % "-"
                                                       % "-"
                                 ).str());
            }
        }
        
        o.writer->close();
    }
}