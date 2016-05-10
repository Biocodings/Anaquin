#include "MetaQuin/m_blat.hpp"
#include "parsers/parser_tsv.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMSen();

/*
 * Implements deno-assembly by a simple implementation of hashing.
 */

struct DAHashImpl : DAsssembly::DenoAssemblyImpl
{
    inline SequinID findSeq(const ContigID &cID) const
    {
        return c2s.count(cID) ? c2s.at(cID) : "";
    }

    std::map<ContigID, SequinID> c2s;
};

MAssembly::Stats MAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_meta;

    DAHashImpl impl;
    MAssembly::Stats stats;

    o.analyze(file);

    switch (o.soft)
    {
        case MetaQuast:
        {
            ParserQuast::parseContigs(Reader(o.contigs), [&](const ParserQuast::ContigData &x, const ParserProgress &)
            {
                for (const auto &c : x.contigs)
                {
                    // Contigs.fasta doesn't have "_"
                    auto t = c;
                    
                    // Eg: contig-1056000000 2818 nucleotides
                    boost::replace_all(t, "_", " ");
                    
                    impl.c2s[t] = x.id;
                }
            });

            ParserQuast::parseGenome(Reader(o.genome), [&](const ParserQuast::GenomeData &x, const ParserProgress &)
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
    
    stats.dnovo = DAsssembly::analyze(file, &impl);
    
    return stats;
}

//static Scripts generateSummary(const FileName &file, const MAssembly::Stats &stats, const MAssembly::Options &o)
//{
//    const auto &r = Standard::instance().r_meta;
//
//    const auto summary = "Summary for input: %1%\n\n"
//                         "   Synthetic: %2%\n"
//                         "   Community: %3%\n\n"
//                         "   Contigs:   %4%\n"
//                         "   Assembled: %5%\n\n"
//                         "   ***\n"
//                         "   *** Reference annotation (Synthetic)\n"
//                         "   ***\n\n"
//                         "   File: %6%\n\n"
//                         "   Synthetic: %7% sequins\n\n"
//                         "   ***\n"
//                         "   ***\n"
//                         "   *** The following statistics are computed on the synthetic community\n"
//                         "   ***\n\n"
//                         "   N20:      %8%\n"
//                         "   N50:      %9%\n"
//                         "   N80:      %10%\n"
//                         "   Min:      %11%\n"
//                         "   Mean:     %12%\n"
//                         "   Max:      %13%\n\n"
//                         "   ***\n"
//                         "   *** The following overlapping statistics are computed as proportion\n"
//                         "   ***\n\n"
//                         "   Match:    %14%\n"
//                         "   Mismatch: %15%\n"
//                         "   Gaps (sequins): %16%\n"
//                         "   Gaps (contigs): %17%\n";
//    
//    return (boost::format(summary) % file
//                                   % stats.blat.n_chrT
//                                   % stats.blat.n_geno
//                                   % stats.blat.aligns.size()
//                                   % stats.blat.countAssembled()
//                                   % o.rChrT
//                                   % r.data().size()
//                                   % stats.N20
//                                   % stats.N50
//                                   % stats.N80
//                                   % stats.min
//                                   % stats.mean
//                                   % stats.max
//                                   % stats.blat.overMatch()
//                                   % stats.blat.overRGaps()
//                                   % stats.blat.overQGaps()
//                                   % stats.blat.overMismatch()).str();
//}

void MAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = MAssembly::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.info("Generating MetaAssembly_summary.stats");
    o.writer->open("MetaAssembly_summary.stats");
    //o.writer->write(generateSummary("MetaAssembly_summary.stats", stats, o));
    o.writer->close();
    
    /*
     * Generating detailed statistics
     */
    
    o.info("Generating MetaAssembly_quins.stats");
    o.writer->open("MetaAssembly_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating limit-of-assembly (LOA)
     */
    
    o.info("Generating MetaAssembly_assembly.R");
    o.writer->open("MetaAssembly_assembly.R");
    o.writer->write(RWriter::createScript("MetaAssembly_quins.stats", PlotMSen()));
    o.writer->close();

    /*
     * Generating detailed statistics for each contig
     */

//    {
//        o.info("Generating MetaAssembly_contigs.stats");
//        o.writer->open("MetaAssembly_contigs.stats");
//        
//        const auto format = "%1%\t%2%";
//        
//        o.writer->write((boost::format(format) % "contig" % "sequin").str());
//        
//        for (const auto &i : stats.blat.aligns)
//        {
//            o.writer->write((boost::format(format) % i.first % i.second->id()).str());
//        }
//        
//        o.writer->close();
//    }

    /*
     * Generating detailed statistics for each sequin
     */

    {
//        o.logInfo("Generating MetaAssembly_quins.stats");
//        o.writer->open("MetaAssembly_quins.stats");
//        
//        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";
//
//        o.writer->write((boost::format(format) % "sequin"
//                                               % "contig"
//                                               % "covered"
//                                               % "match"
//                                               % "mismatch"
//                                               % "tgap"
//                                               % "qgap").str());
//
//        for (const auto &i : stats.blat.metas)
//        {
//            const auto &align = i.second;
//            
//            o.writer->write((boost::format(format) % align->seq->id
//                                                   % align->contigs.size()
//                                                   % align->covered
//                                                   % align->oMatch
//                                                   % align->oMismatch
//                                                   % align->oRGaps
//                                                   % align->oQGaps).str());
//        }
//        
//        o.writer->close();
    }
}