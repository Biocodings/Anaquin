#include "MetaQuin/m_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMReads();

MAlign::Stats MAlign::analyze(const FileName &file, const MAlign::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MAlign::Stats stats;
    stats.hist = r.hist();

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &)
    {
        if (align.mapped)
        {
            const auto m = r.match(align.cID);
            
            if (m)
            {
                stats.n_chrT++;
                stats.hist.at(m->id)++;
            }
            else
            {
                stats.n_geno++;
            }
        }
        else
        {
            stats.n_unmap++;
        }
    });
   
    /*
     * Build a linear model for input concentration against number of reads (on the logarithm scale)
     */
    
    for (auto &i : stats.hist)
    {
        if (i.second)
        {
            stats.add(i.first, r.match(i.first)->concent(), i.second);
        }
    }

    stats.limit = r.absolute(stats.hist);

    return stats;
}

static void writeQuins(const FileName &file, const MAlign::Stats &stats, const MAlign::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);

    o.writer->open(file);
    o.writer->write((boost::format(format) % "seq" % "input" % "reads").str());
    
    const auto total = sum(stats.hist);
    
    for (const auto &i : stats.hist)
    {
        const auto l = r.match(i.first)->l;
        assert(l.length() > 1);
        
        // Input concentration (attomol/ul)
        const auto expected = r.match(i.first)->concent();
        
        // Measured FPKM
        const auto measured = ((double)i.second * pow(10, 9)) / (total * l.length());

        o.writer->write((boost::format(format) % i.first % expected % measured).str());
    }

    o.writer->close();
}

void MAlign::report(const FileName &file, const MAlign::Options &o)
{
    const auto stats = MAlign::analyze(file, o);

    /*
     * Generating MetaAlign_summary.stats
     */

    o.generate("MetaAlign_summary.stats");
    o.writer->open("MetaAlign_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rGeno,
                                                file,
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();

    /*
     * Generating MetaAlign_quins.stats
     */
    
    writeQuins("MetaAlign_quins.stats", stats, o);
    
    /*
     * Generating MetaAlign_reads.R
     */
    
    o.generate("MetaAlign_reads.R");
    o.writer->open("MetaAlign_reads.R");
    o.writer->write(RWriter::createScript("MetaAlign_quins.stats", PlotMReads()));
    o.writer->close();
}