#include "data/types.hpp"
#include "parsers/parser_tsv.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMAssembly();

MAssembly::Stats MAssembly::analyze(const std::vector<FileName> &files, const Options &o)
{
    // Eg: Contigs.fasta
    const auto fasta = files[0];

    // Eg: align.psl and alignments_Contigs.tsv
    const auto align = files[1];
    
    const auto &r = Standard::instance().r_meta;

    MAssembly::Stats stats;

    o.analyze(fasta);

    /*
     * Constructing mapping involving contigs and sequins. However, how this is constructed is tool specific.
     */
    
    switch (o.align)
    {
        case Blat:
        {
            const auto x = MBlat::analyze(align);
            
            for (auto &i : x.aligns)
            {
                stats.c2s[i.first] = i.second->id();
            }

            for (auto &i : x.metas)
            {
                for (auto &j : i.second->contigs)
                {
                    stats.s2c[i.first].push_back(j.id);
                    stats.add(i.first, 10, 20); // TODO: Fix this
                }
            }
            
            /*
             * TODO: Calculating the sensitivity for each individual sequin
             */
            
            //for (auto &i : x.s2c)
            {
                //stats.add(i.first, 10, 20);
            }

            break;
        }

        case MetaQuast:
        {
            ParserQuast::parseAlign(Reader(align), [&](const ParserQuast::ContigData &x,
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

            // Eg: genome_info.txt
            const auto genome = files[2];
            
            ParserQuast::parseGenome(Reader(genome), [&](const ParserQuast::GenomeData &x,
                                                         const ParserProgress &)
            {
                const auto match = r.match(x.id);

                if (match)
                {
                    // Build a linear model between input concentration and sensitivity
                    stats.add(match->id, match->concent(), static_cast<Proportion>(x.covered) / x.total);
                }
            });

            break;
        }
    }
    
    assert(!stats.c2s.empty());
    assert(!stats.s2c.empty());

    stats.soft  = o.soft;
    stats.align = o.align;
    
    // Calculate statistics such as N50 and proportion asssembled
    stats.dnovo = DAsssembly::analyze(fasta, &stats);

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
                         "   *** The following statistics are computed on the synthetic community\n"
                         "   ***\n\n"
                         "   N20:  %8%\n"
                         "   N50:  %9%\n"
                         "   N80:  %10%\n"
                         "   Min:  %11%\n"
                         "   Mean: %12%\n"
                         "   Max:  %13%\n\n"
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

static Scripts writeContigs(const MAssembly::Stats &stats, const MAssembly::Options &o)
{
    const auto &r = Standard::instance().r_meta;

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";

    std::stringstream ss;
    ss << ((boost::format(format) % "seq"
                                  % "input"
                                  % "contig"
                                  % "match"
                                  % "mismatch")) << std::endl;

    for (const auto &seq : r.data())
    {
        if (stats.s2c.count(seq.first))
        {
            for (const auto &c : stats.s2c.at(seq.first))
            {
                switch (o.align)
                {
                    case MAssembly::Blat:
                    {
                        throw "Not Implementd";
                        break;
                    }
                        
                    case MAssembly::MetaQuast:
                    {
                        /*
                         * The alignment input: "genome_info.txt" combines the sensitivity for all contigs aligned
                         * to the sequin. Thus, it's not possible to give the information at the contig level.
                         */
                        
                        ss << ((boost::format(format) % seq.first
                                                      % seq.second.concent()
                                                      % c
                                                      % "-"
                                                      % "-").str());
                        ss << std::endl;

                        break;
                    }
                }
                
            }
        }
    }

    return ss.str();
}

void MAssembly::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = MAssembly::analyze(files, o);

    o.info("Generating statistics");
    
    /*
     * Generating MetaAssembly_summary.stats
     */

    o.info("Generating MetaAssembly_summary.stats");
    o.writer->open("MetaAssembly_summary.stats");
    o.writer->write(generateSummary("MetaAssembly_summary.stats", stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_quins.stats
     */
    
    o.info("Generating MetaAssembly_quins.stats");
    o.writer->open("MetaAssembly_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating MetaAssembly_queries.stats
     */
    
    o.info("Generating MetaAssembly_queries.stats");
    o.writer->open("MetaAssembly_queries.stats");
    o.writer->write(writeContigs(stats, o));
    o.writer->close();

    /*
     * Generating MetaAssembly_assembly.R
     */

    o.info("Generating MetaAssembly_assembly.R");
    o.writer->open("MetaAssembly_assembly.R");
    o.writer->write(RWriter::createScript("MetaAssembly_quins.stats", PlotMAssembly()));
    o.writer->close();
}