#include "data/data.hpp"
#include "parsers/parser_tsv.hpp"
#include "MetaQuin/m_assembly.hpp"

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
            
            // For each sequin...
            for (auto &s : x.metas)
            {
                // Sequin name
                const auto &id = s.first;
                
                // Length of the sequin
                const auto l = s.second->seq.l;
                
                // Required for detecting overlapping
                DInter i(s.first, Locus(0, l.length()));
                
                /*
                 * For each contig aligned to the sequin, we calculate where exactly it aligns.
                 */
                
                for (auto &j : s.second->contigs)
                {
                    stats.s2c[s.first].push_back(j.id);
                    i.add(Locus(j.l.start, j.l.end - 1));
                }
                
                // Statistics after accounting for all mapping contigs
                const auto si = i.stats();

                assert(si.nonZeros + si.zeros == si.length);
                
                stats.z  += si.zeros;
                stats.nz += si.nonZeros;
                
                // Number of bases with coverage for the squin
                stats.s2nz[id] = si.nonZeros;
                
                stats.add(id, r.input1(id, o.mix), si.covered());
            }
            
            break;
        }

        default : { throw "Not Implemented"; }
    }
    
    A_ASSERT(!stats.c2s.empty());
    A_ASSERT(!stats.s2c.empty());
    
    stats.dn = DAsssembly::analyze(fasta, &stats);
    
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
        else
        {
            o.warn(seq + " not assembled");
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
                         "       Match:       %13%\n"
                         "       Mismatch:    %14%\n"
                         "       Covered:     %15% bases\n"
                         "       Not covered: %16% bases\n"
                         "       Sensitivity: %17%\n";

    const auto &dn = stats.dn;
    
    // Overall sensitivity
    const auto sn = (Proportion) stats.nz / (stats.z + stats.nz);
    
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
                                   % stats.z
                                   % stats.nz
                                   % sn).str();
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
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    ss << (boost::format(format) % "Name"
                                 % "Length"
                                 % "Input"
                                 % "Covered"
                                 % "Sn").str() << std::endl;
    
    const auto r1 = Standard::instance().r_meta.regs1();
    
    for (const auto &i : stats)
    {
        const auto &id = i.first;
        ss << ((boost::format(format) % id
                                      % r1.at(id).length()
                                      % i.second.x
                                      % (stats.s2nz.count(id) ? stats.s2nz.at(id) : 0)
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
    
    o.info("Generating MetaAssembly_sequins.tsv");
    o.writer->open("MetaAssembly_sequins.tsv");
    o.writer->write(MAssembly::generateQuins(stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_queries.csv
     */
    
    o.info("Generating MetaAssembly_queries.tsv");
    o.writer->open("MetaAssembly_queries.tsv");
    o.writer->write(writeContigs(stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAssembly_assembly.R
     */

    o.info("Generating MetaAssembly_assembly.R");
    o.writer->open("MetaAssembly_assembly.R");
    o.writer->write(RWriter::createLogistic("MetaAssembly_sequins.tsv",
                                            "Assembly Detection",
                                            "Input Concentration (log2)",
                                            "Sensitivity",
                                            "Input",
                                            "Sn",
                                            true));
    o.writer->close();
}
