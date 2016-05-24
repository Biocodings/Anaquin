#include "data/types.hpp"
#include "parsers/parser_tsv.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMAssembly();

MAssembly::Stats MAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_meta;

    MAssembly::Stats stats;

    o.analyze(file);

    switch (o.soft)
    {
        case MetaQuast:
        {
            ParserQuast::parseContigs(Reader(o.contigs), [&](const ParserQuast::ContigData &x,
                                                             const ParserProgress &)
            {
                for (const auto &c : x.contigs)
                {
                    // Contigs.fasta doesn't have "_"
                    auto t = c;
                    
                    // Eg: contig-1056000000 2818 nucleotides
                    boost::replace_all(t, "_", " ");
                    
                    stats.c2s[t] = x.id;
                    stats.s2c[x.id].push_back(t);
                }
            });

            ParserQuast::parseGenome(Reader(o.genome), [&](const ParserQuast::GenomeData &x,
                                                           const ParserProgress &)
            {
                const auto match = r.match(x.id);

                if (match)
                {
                    stats.add(match->id, match->concent(), static_cast<Proportion>(x.covered) / x.total);
                }
            });

            break;
        }
    }
    
    stats.dnovo = DAsssembly::analyze(file, &stats);

    return stats;
}

static Scripts generateSummary(const FileName &file, const MAssembly::Stats &stats, const MAssembly::Options &o)
{
    const auto &r = Standard::instance().r_meta;

    const auto summary = "Summary for input: %1%\n\n"
                         "   Synthetic: %2%\n"
                         "   Community: %3%\n\n"
                         "   Contigs:   %4%\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %6%\n\n"
                         "   Synthetic: %7% sequins\n"
                         "   Contigs: %5%\n\n"
                         "   ***\n"
                         "   ***\n"
                         "   *** The following statistics are computed on the synthetic community\n"
                         "   ***\n\n"
                         "   N20:      %8%\n"
                         "   N50:      %9%\n"
                         "   N80:      %10%\n"
                         "   Min:      %11%\n"
                         "   Mean:     %12%\n"
                         "   Max:      %13%\n\n"
                         "   ***\n"
                         "   *** The following overlapping statistics are computed as proportion\n"
                         "   ***\n\n"
                         "   Match:    %14%\n"
                         "   Mismatch: %15%\n"
                         "   Gaps (sequins): %16%\n"
                         "   Gaps (contigs): %17%\n";
    
    const auto &dn = stats.dnovo;
    
    return (boost::format(summary) % file
                                   % dn.n_chrT
                                   % dn.n_geno
                                   % (dn.n_chrT + dn.n_geno)
                                   % dn.n_chrT
                                   % o.rChrT
                                   % r.data().size()
                                   % dn.N20
                                   % dn.N50
                                   % dn.N80
                                   % dn.min
                                   % dn.mean
                                   % dn.max
                                   % "-" //stats.blat.overMatch()
                                   % "-" //stats.blat.overRGaps()
                                   % "-" //stats.blat.overQGaps()
                                   % "-" //stats.blat.overMismatch())
            ).str();
}

static Scripts generateQuins(const MAssembly::Stats &stats)
{
    const auto &r = Standard::instance().r_meta;

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

    std::stringstream ss;
    ss << ((boost::format(format) % "seq"
                                  % "input"
                                  % "contig"
                                  % "covered"
                                  % "match"
                                  % "mismatch"
                                  % "tgap"
                                  % "qgap")) << std::endl;

    // Statistics for sensitivity
    const auto sst = stats.data(false);

    for (const auto &seq : r.data())
    {
        if (stats.s2c.count(seq.first))
        {
            for (const auto &c : stats.s2c.at(seq.first))
            {
                ss << ((boost::format(format) % seq.first
                                              % seq.second.concent()
                                              % c
                                              % sst.id2y.at(seq.first)
                                              % "-"
                                              % "-"
                                              % "-"
                                              % "-").str());
            }
        }
        else
        {
            ss << ((boost::format(format) % seq.first
                                          % seq.second.concent()
                                          % "-"
                                          % "-"
                                          % "-"
                                          % "-"
                                          % "-"
                                          % "-").str());
        }
    }

    return ss.str();
}

static Scripts generateContigs(const MAssembly::Stats &stats)
{
    const auto format = "%1%\t%2%";
    
    std::stringstream ss;
    ss << (boost::format(format) % "contig" % "seq").str() << std::endl;
    
    for (const auto &i : stats.c2s)
    {
        ss << (boost::format(format) % i.first % i.second).str() << std::endl;
    }
    
    return ss.str();
}

void MAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = MAssembly::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.info("Generating MetaAssembly_summary.stats");
    o.writer->open("MetaAssembly_summary.stats");
    o.writer->write(generateSummary("MetaAssembly_summary.stats", stats, o));
    o.writer->close();
    
    /*
     * Generating detailed statistics
     */
    
    o.info("Generating MetaAssembly_quins.stats");
    o.writer->open("MetaAssembly_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating MetaAssembly_contigs.stats
     */

    o.info("Generating MetaAssembly_contigs.stats");
    o.writer->open("MetaAssembly_contigs.stats");
    o.writer->write(generateContigs(stats));
    o.writer->close();

    /*
     * Generating MetaAssembly_quins.stats
     */

    o.info("Generating MetaAssembly_quins.stats_");
    o.writer->open("MetaAssembly_quins.stats_");
    o.writer->write(generateQuins(stats));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_assembly.R
     */

    o.info("Generating MetaAssembly_assembly.R");
    o.writer->open("MetaAssembly_assembly.R");
    o.writer->write(RWriter::createScript("MetaAssembly_quins.stats", PlotMAssembly()));
    o.writer->close();
}