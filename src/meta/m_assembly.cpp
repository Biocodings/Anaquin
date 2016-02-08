#include "meta/m_blat.hpp"
#include "meta/m_assembly.hpp"
#include "parsers/parser_tsv.hpp"

using namespace Anaquin;

MAssembly::Stats MAssembly::analyze(const FileName &file, const Options &o)
{
    MAssembly::Stats stats;

    assert(!o.psl.empty());
    assert(o.tool != MetaAssembler::RayMeta || !o.contigs.empty());

    /*
     * Generate statistics for the alignment
     */
    
    o.analyze(o.psl);
    
    // Analyse the blat alignment file
    const auto t = MBlat::analyze(o.psl);
 
    o.info("Found: " + std::to_string(t.aligns.size()) + " in " + o.psl);
    
    /*
     * Generate statistics for the assembler
     */

    o.analyze(file);

    switch (o.tool)
    {
        case Velvet:  { stats = Velvet::analyze<MAssembly::Stats, Contig>(file, &t);             break; }
        case RayMeta: { stats = RayMeta::analyze<MAssembly::Stats, Contig>(file, o.contigs, &t); break; }
    }
    
    stats.blat = t;

    return stats;
}

void MAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = MAssembly::analyze(file, o);

    /*
     * Generating summary statistics
     */

    {
        o.logInfo("Generating summary statistics");
        o.writer->open("MetaAssembly_summary.stats");
        
        const auto summary = "Summary for input: %1%\n\n"
                             "   Community: %2%\n"
                             "   Synthetic: %3%\n\n"
                             "   Contigs:   %4%\n"
                             "   Assembled: %5%\n"
                             "   Reference: %6%\n\n"
                             "   ***\n"
                             "   *** The following statistics are computed on the synthetic community\n"
                             "   ***\n\n"
                             "   N20:      %7%\n"
                             "   N50:      %8%\n"
                             "   N80:      %9%\n"
                             "   Min:      %10%\n"
                             "   Mean:     %11%\n"
                             "   Max:      %12%\n\n"
                             "   ***\n"
                             "   *** The following overlapping statistics are computed as proportion\n"
                             "   ***\n\n"
                             "   Match:    %13%\n"
                             "   Mismatch: %14%\n"
                             "   Gaps (sequins): %15%\n"
                             "   Gaps (contigs): %16%\n";

        o.writer->write((boost::format(summary) % file
                                                % stats.blat.n_endo
                                                % stats.blat.n_chrT
                                                % stats.blat.aligns.size()
                                                % stats.blat.countAssembled()
                                                % stats.blat.metas.size()
                                                % stats.N20
                                                % stats.N50
                                                % stats.N80
                                                % stats.min
                                                % stats.mean
                                                % stats.max
                                                % stats.blat.overMatch()
                                                % stats.blat.overRGaps()
                                                % stats.blat.overQGaps()
                                                % stats.blat.overMismatch()).str());
        o.writer->close();
    }

    /*
     * Generating detailed statistics for each contig
     */

    {
        o.writer->open("MetaAssembly_contigs.stats");
        
        const std::string format = "%1%\t%2%";
        
        o.writer->write((boost::format(format) % "ID" % "Sequin ID").str());
        
        for (const auto &i : stats.blat.aligns)
        {
            o.writer->write((boost::format(format) % i.first % i.second->id()).str());
        }

        o.writer->close();
    }
    
    /*
     * Generating detailed statistics for each sequin
     */

    {
        o.logInfo("Generating sequins statistics");
        o.writer->open("MetaAssembly_quins.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";

        o.writer->write((boost::format(format) % "ID"
                                               % "Contigs"
                                               % "Covered"
                                               % "Match"
                                               % "Mismatch"
                                               % "TGaps"
                                               % "QGaps").str());

        for (const auto &i : stats.blat.metas)
        {
            const auto &align = i.second;
            
            o.writer->write((boost::format(format) % align->seq->id
                                                   % align->contigs.size()
                                                   % align->covered
                                                   % align->oMatch
                                                   % align->oMismatch
                                                   % align->oRGaps
                                                   % align->oQGaps).str());
        }
        
        o.writer->close();
    }
}