#include "data/standard.hpp"
#include "MetaQuin/m_abund.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMReads();

MAbund::Stats MAbund::analyze(const FileName &file, const MAbund::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MAbund::Stats stats;
    stats.hist = r.hist();

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &)
    {
        if (align.mapped)
        {
            const auto m = r.match(align.cID);
            
            if (m)
            {
                stats.n_syn++;
                stats.hist.at(m->id)++;
            }
            else
            {
                stats.n_gen++;
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

//    stats.limit = r.detectLimit(stats.hist);

    return stats;
}

static void writeQuins(const FileName &file, const MAbund::Stats &stats, const MAbund::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);

    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID" % "input" % "reads").str());
    
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

void MAbund::report(const FileName &file, const MAbund::Options &o)
{
    const auto stats = MAbund::analyze(file, o);

    /*
     * Generating MetaAbund_summary.stats
     */

    o.generate("MetaAbund_summary.stats");
    o.writer->open("MetaAbund_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rAnnot,
                                                o.rAnnot,
                                                file,
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();

    /*
     * Generating MetaAbund_quins.stats
     */
    
    writeQuins("MetaAbund_sequins.csv", stats, o);
}