#include "data/standard.hpp"
#include "MetaQuin/m_abund.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

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
                stats.nSyn++;
                stats.hist.at(m->id)++;
                
                if (stats.limit.id.empty() || m->concent() < stats.limit.abund)
                {
                    stats.limit.id = m->id;
                    stats.limit.abund = m->concent();
                }
            }
            else
            {
                stats.nGen++;
            }
        }
        else
        {
            stats.nNA++;
        }
    });

    for (auto &i : stats.hist)
    {
        if (i.second)
        {
            stats.add(i.first, r.match(i.first)->concent(), i.second);
        }
    }
    
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
        
        // Input concentration (attomol/ul)
        const auto expected = r.match(i.first)->concent();
        
        // Measured FPKM
        const auto measured = ((double)i.second * pow(10, 9)) / (total * l.length());
        
        o.writer->write((boost::format(format) % i.first % expected % measured).str());
    }
    
    o.writer->close();
}

Scripts MAbund::generateRLinear(const FileName &src, const Stats &stats, const Options &o)
{
    return RWriter::createRLinear(src,
                                  o.work,
                                  "Allele Frequency",
                                  "Expected allele frequency (log2)",
                                  "Measured allele frequency (log2)",
                                  "log2(data$ExpFreq)",
                                  "log2(data$ObsFreq)",
                                  "input",
                                  true);
}

static void writeRLinear(const FileName &file, const FileName &src, const MAbund::Stats &stats, const MAbund::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(MAbund::generateRLinear(src, stats, o));
    o.writer->close();
}

static Scripts generateSummary(const FileName &src, const MAbund::Stats &stats, const MAbund::Options &o)
{
    // Defined in resources.cpp
    extern FileName MixRef();

    const auto &r = Standard::instance().r_meta;
    const auto ls = stats.linear();
    
    const auto format = "-------MetaAbund Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference MetaQuin Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------Sequin Counts\n\n"
                        "       Synthetic: %4%\n"
                        "       Detection Sensitivity: %5% (attomol/ul) (%6%)\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %7%\n"
                        "       Correlation: %8%\n"
                        "       R2:          %9%\n"
                        "       F-statistic: %10%\n"
                        "       P-value:     %11%\n"
                        "       SSM:         %12%, DF: %13%\n"
                        "       SSE:         %14%, DF: %15%\n"
                        "       SST:         %16%, DF: %17%\n";
    
    return (boost::format(format) % src               // 1
                                  % r.countSeqs()     // 2
                                  % MixRef()          // 3
                                  % stats.size()      // 4
                                  % stats.limit.abund // 5
                                  % stats.limit.id    // 6
                                  % ls.m              // 7
                                  % ls.r              // 8
                                  % ls.R2             // 9
                                  % ls.F              // 10
                                  % ls.p              // 11
                                  % ls.SSM            // 12
                                  % ls.SSM_D          // 13
                                  % ls.SSE            // 14
                                  % ls.SSE_D          // 15
                                  % ls.SST            // 16
                                  % ls.SST_D          // 17
            ).str();
}

void MAbund::report(const FileName &file, const MAbund::Options &o)
{
    const auto stats = MAbund::analyze(file, o);
    
    /*
     * Generating MetaAbund_summary.stats
     */
    
    o.generate("MetaAbund_summary.stats");
    o.writer->open("MetaAbund_summary.stats");
    o.writer->write(generateSummary(file, stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAbund_sequins.csv
     */
    
    writeQuins("MetaAbund_sequins.csv", stats, o);
    
    /*
     * Generating MetaAbund_linear.R
     */
    
    writeRLinear("MetaAbund_linear.R", "MetaAbund_sequins.csv", stats, o);
}
