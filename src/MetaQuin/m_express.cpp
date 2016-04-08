#include "MetaQuin/m_express.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAbundAbund();

MExpress::Stats MExpress::analyze(const FileName &file, const MExpress::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MExpress::Stats stats;

    // Initialize the sequins
    stats.hist = r.hist();
    
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

    switch (o.soft)
    {
        case Software::Velvet:  { stats.assembly = Velvet::analyze<MAssembly::Stats, Contig>(file, &t);             break; }
        case Software::RayMeta: { stats.assembly = RayMeta::analyze<MAssembly::Stats, Contig>(file, o.contigs, &t); break; }
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
    stats.n_endo = stats.assembly.n - stats.n_chrT;

    o.info("Analyzing the alignments");

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

        if (stats.limit.id.empty() || align->seq->abund(Mix_1, false) < stats.limit.abund)
        {
            stats.limit.id     = align->seq->id;
            stats.limit.abund  = align->seq->abund(Mix_1, false);
            stats.limit.counts = align->contigs.size();
        }
        
        const auto p = MExpress::calculate(stats, stats.blat, stats.assembly, align->seq->id, *meta.second, o, o.coverage);
        
        if (p.x && p.y)
        {
            stats.add(align->seq->id, p.x, p.y);
        }
    }

    stats.limit = r.absolute(stats.hist);

    return stats;
}

static void generateContigs(const FileName &file, const MExpress::Stats &stats, const MExpress::Options &o)
{
    o.info("Generating " + file);
    o.writer->open(file);
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.writer->write((boost::format(format) % "contigID"
                                           % "seqID"
                                           % "length"
                                           % "coverage"
                                           % "normalized").str());
    
    for (const auto &i : stats.blat.aligns)
    {
        if (stats.assembly.contigs.count(i.first))
        {
            const auto &contig = stats.assembly.contigs.at(i.first);
            
            o.writer->write((boost::format(format) % i.first
                                                   % i.second->id()
                                                   % contig.k_len
                                                   % contig.k_cov
                                                   % contig.normalized()).str());
        }
        else
        {
            o.writer->write((boost::format(format) % i.first
                                                   % i.second->id()
                                                   % "-"
                                                   % "-"
                                                   % "-").str());
        }
    }
    
    o.writer->close();
}

void MExpress::report(const FileName &file, const MExpress::Options &o)
{
    const auto stats = MExpress::analyze(file, o);

    /*
     * 1. Generating summary statistics
     */
    
    o.info("Generating MetaAbund_summary.stats");
    o.writer->open("MetaAbund_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                file,
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();
    
    /*
     * 2. Generating CSV for all sequins
     */
    
    o.info("Generating MetaAbund_quins.csv");
    o.writer->open("MetaAbund_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * 3. Generating for AbundAbund
     */
    
    o.info("Generating MetaAbund_abundAbund.R");
    o.writer->open("MetaAbund_abundAbund.R");
    o.writer->write(RWriter::createScript("MetaAbund_quins.csv", PlotVAbundAbund()));
    o.writer->close();
    
    /*
     * 4. Generating detailed statistics for the contigs
     */
    
    generateContigs("MetaAbund_contigs.stats", stats, o);
}