#include "data/data.hpp"
#include "parsers/parser_tsv.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

MAssembly::Stats MAssembly::analyze(const std::vector<FileName> &files, const Options &o)
{
    // Eg: Contigs.fasta
    const auto fasta = files[0];
    
    // Eg: align.psl or alignments_Contigs.tsv
    const auto align = files[1];
    
    const auto &r = Standard::instance().r_meta;
    
    MAssembly::Stats stats;
    
    o.analyze(fasta);
    
    switch (o.format)
    {
        case Format::Blat:
        {
            const auto x = MBlat::analyze(align);
            
            /*
             * Building mapping for contigs
             */
            
            stats.c2s = x.c2s;
            stats.c2l = x.c2l;
            stats.c2a = x.c2a;
            
            A_ASSERT(!stats.c2s.empty());
            A_ASSERT(!stats.c2l.empty());
            A_ASSERT(!stats.c2a.empty());
            
            /*
             * Building mapping for sequins
             */
            
            for (auto &s : x.metas)
            {
                // Length of the sequin
                const auto l = s.second->seq.l;
                
                // Required for detecting overlapping
                DInter i(s.first, Locus(0, l.length()));
                
                for (auto &j : s.second->contigs)
                {
                    stats.s2c[s.first].push_back(j.id);
                    i.add(Locus(j.l.start, j.l.end - 1));
                }

                stats.add(s.first, r.input1(s.first, o.mix), i.stats().covered());
            }
            
            break;
        }

        default : { throw "Not Implemented"; }
    }
    
    A_ASSERT(!stats.c2s.empty());
    A_ASSERT(!stats.s2c.empty());
    
    // Calculate statistics such as N50 and proportion asssembled
    stats.dnovo = DAsssembly::analyze(fasta, &stats);
    
    /*
     * Calculating the proportion being assembled
     */
    
    for (const auto &seq : r.seqsL1())
    {
        if (stats.s2c.count(seq))
        {
            for (const auto &c : stats.s2c.at(seq))
            {
                switch (o.format)
                {
                    case Format::Blat:
                    {
                        stats.match += stats.c2a.at(c);
                        stats.mismatch += (stats.c2l.at(c) - stats.c2a.at(c));
                        break;
                    }
                }
            }
        }
    }
    
    return stats;
}

static Scripts generateSummary(const FileName &src, const MAssembly::Stats &stats, const MAssembly::Options &o)
{
    extern FileName BedRef();
    
    const auto summary = "-------MetaAssembly Output\n\n"
                         "       Summary for input: %1%\n\n"
                         "       Synthetic: %2% contigs\n"
                         "       Genome:    %3% contigs\n"
                         "       Total:     %4% contigs\n\n"
                         "-------Reference MetaQuin Annotation\n\n"
                         "       File: %5%\n"
                         "       Synthetic: %6% sequins\n\n"
                         "-------Assembly Statistics\n\n"
                         "       ***\n"
                         "       *** The following statistics are computed on the synthetic community\n"
                         "       ***\n\n"
                         "       N20:  %7%\n"
                         "       N50:  %8%\n"
                         "       N80:  %9%\n"
                         "       Min:  %10%\n"
                         "       Mean: %11%\n"
                         "       Max:  %12%\n\n"
                         "       ***\n"
                         "       *** The following overlapping statistics are computed as proportion\n"
                         "       ***\n\n"
                         "       Match:       %13%\n"
                         "       Mismatch:    %14%\n"
                         "       Sensitivity: %15%\n";

    const auto &dn = stats.dnovo;
    
    return (boost::format(summary) % src
                                   % dn.nSeqs
                                   % dn.nEndo
                                   % (dn.nSeqs + dn.nEndo)
                                   % BedRef()
                                   % Standard::instance().r_meta.seqsL1().size()
                                   % dn.N20
                                   % dn.N50
                                   % dn.N80
                                   % dn.min
                                   % dn.mean
                                   % dn.max
                                   % stats.match
                                   % stats.mismatch
                                   % stats.covered()).str();
}

static Scripts writeContigs(const MAssembly::Stats &stats, const MAssembly::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    std::stringstream ss;
    ss << ((boost::format(format) % "Name"
                                  % "Input"
                                  % "Contig"
                                  % "Match"
                                  % "Mismatch")) << std::endl;
    
    for (const auto &seq : r.seqsL1())
    {
        if (stats.s2c.count(seq))
        {
            for (const auto &c : stats.s2c.at(seq))
            {
                switch (o.format)
                {
                    case MAssembly::Format::Blat:
                    {
                        const auto total = stats.c2l.at(c);
                        const auto align = stats.c2a.at(c);
                        
                        assert(total >= align);
                        
                        ss << ((boost::format(format) % seq
                                                      % r.input1(seq, o.mix)
                                                      % c
                                                      % align
                                                      % (total - align)).str()) << std::endl;
                        break;
                    }
                }
            }
        }
    }
    
    return ss.str();
}

Scripts MAssembly::generateQuins(const Stats &stats, const Options &o)
{
    std::stringstream ss;
    
    const auto format = "%1%\t%2%\t%3%\t%4%";
    ss << (boost::format(format) % "Name" % "Length" % "Input" % "Sn").str() << std::endl;
    
    const auto r1 = Standard::instance().r_meta.regs1();
    
    for (const auto &i : stats)
    {
        ss << ((boost::format(format) % i.first
                                      % r1.at(i.first).length()
                                      % i.second.x
                                      % i.second.y).str()) << std::endl;
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
    o.writer->write(generateSummary(files[0], stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_sequins.csv
     */
    
    o.info("Generating MetaAssembly_sequins.csv");
    o.writer->open("MetaAssembly_sequins.csv");
    o.writer->write(MAssembly::generateQuins(stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_queries.csv
     */
    
    o.info("Generating MetaAssembly_queries.csv");
    o.writer->open("MetaAssembly_queries.csv");
    o.writer->write(writeContigs(stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_assembly.R
     */

    o.info("Generating MetaAssembly_assembly.R");
    o.writer->open("MetaAssembly_assembly.R");
    o.writer->write(RWriter::createLogistic("MetaAssembly_sequins.csv",
                                            "Assembly Detection",
                                            "Input Concentration (log2)",
                                            "Sensitivity",
                                            "Input",
                                            "Sn",
                                            true));
    o.writer->close();
}
