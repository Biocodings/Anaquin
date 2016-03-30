#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats MAbundance::analyze(const FileName &file, const MAbundance::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
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
        case Velvet:  { stats.chrT->assembly = Velvet::analyze<MAssembly::Stats, Contig>(file, &t);             break; }
        case RayMeta: { stats.chrT->assembly = RayMeta::analyze<MAssembly::Stats, Contig>(file, o.contigs, &t); break; }
    }

    stats.chrT->blat = t;

    if (!stats.chrT->assembly.n)
    {
        throw std::runtime_error("No contig detected in the input file. Please check and try again.");
    }
    else if (stats.chrT->assembly.contigs.empty())
    {
        throw std::runtime_error("No contig aligned in the input file. Please check and try again.");
    }

    stats.chrT->n_chrT = stats.chrT->assembly.contigs.size();
    stats.chrT->n_endo = stats.chrT->assembly.n - stats.chrT->n_chrT;

    o.info("Analyzing the alignments");

    for (auto &meta : stats.chrT->blat.metas)
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

        if (stats.chrT->limit.id.empty() || align->seq->abund(Mix_1, false) < stats.chrT->limit.abund)
        {
            stats.chrT->limit.id     = align->seq->id;
            stats.chrT->limit.abund  = align->seq->abund(Mix_1, false);
            stats.chrT->limit.counts = align->contigs.size();
        }
        
        const auto p = MAbundance::calculate(stats, stats.chrT->blat, stats.chrT->assembly, align->seq->id, *meta.second, o, o.coverage);
        
        if (p.x && p.y)
        {
            stats.chrT->add(align->seq->id, p.x, p.y);
        }
    }

    stats.chrT->absolute = r.absolute(stats.chrT->h);

    return stats;
}

void MAbundance::report(const FileName &file, const MAbundance::Options &o)
{
    const auto stats = MAbundance::analyze(file, o);

    o.info("Generating summary statistics");
    //AnalyzeReporter::linear("MetaAbund_summary.stats", file, stats, "contigs", o.writer, "sequins");
 
    o.info("Generating scatter plot");
    //AnalyzeReporter::scatter(stats,
                    //         "Expected abundance vs Measured coverage",
                      //       "MetaAbundance",
                        //     "Expected abudnance (attomol/ul)",
                          //   "Measured coverage (k-mer)",
                            // "Expected abdunance (log2 attomol/ul)",
                            // "Measured coverage (log2 k-mer)",
                            // o.writer);
    
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
        
        for (const auto &i : stats.chrT->blat.aligns)
        {
            if (stats.chrT->assembly.contigs.count(i.first))
            {
                const auto &contig = stats.chrT->assembly.contigs.at(i.first);
                
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