#include "trans/t_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

static TAlign::Stats init()
{
    const auto &r = Standard::instance().r_trans;
    
    TAlign::Stats stats;

    // Initalize the distributions
    stats.pb.h = stats.pe.h = stats.pi.h = r.geneHist();

    for (const auto &i : stats.pe.h)
    {
        stats.sb[i.first];
        stats.se[i.first];
        stats.si[i.first];
    }
    
    stats.eInters = r.exonInters();
    stats.iInters = r.intronInters();

    return stats;
}

typedef std::map<GenomeID, Base> FPStats;

static const Interval * matchExon(const Alignment &align, TAlign::Stats &stats, FPStats &fps)
{
    TransRef::ExonInterval * match = nullptr;

    // Can we find an exon that contains the alignment?
    if (classify(stats.pe.m, align, [&](const Alignment &)
    {
        return (match = stats.eInters.contains(align.l));
    }))
    {
        // Update the statistics for the sequin
        classifyTP(stats.se.at(match->gID), align);

        stats.pe.h.at(match->gID)++;
    }

    // Can we find an exon that overlaps the alignment?
    else if ((match = stats.eInters.overlap(align.l)))
    {
        // Update the statistics for the sequin
        classifyFP(stats.se.at(match->gID), align);
    }
    
    if (match)
    {
        // Anything that fails to being mapped is counted as FP
        fps.at(match->gID) += match->map(align.l);
    }
    
    return match;
}

static const Interval * matchIntron(const Alignment &align, TAlign::Stats &stats)
{
    TransRef::IntronInterval * match = nullptr;
    
    if (classify(stats.pi.m, align, [&](const Alignment &)
    {
        return (match = stats.iInters.contains(align.l));
    }))
    {
        // Update the statistics for the sequin
        classifyTP(stats.si.at(match->gID), align);

        stats.pi.h.at(match->gID)++;
    }
    else if ((match = stats.iInters.overlap(align.l)))
    {
        // Update the statistics for the sequin
        classifyFP(stats.si.at(match->gID), align);
    }
    
    return match;
}

TAlign::Stats TAlign::stats(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    TAlign::Stats stats = init();

    o.analyze(file);

    /*
     * Loop through each of the alignments and calculate various metrics at the exon, intron and
     * base level.
     *
     * Exon level is defined as the performance per exon. An alignment that is not mapped entirely
     * within an exon is considered as a FP. The intron level is similar.
     *
     * Base level is defined as the performance per nucleotide. A partial mapped read will have
     * FP and TP.
     *
     * If we can find an exact match, this is obviously a TP. Otherwise, if we
     */

    FPStats fps = stats.pb.h;
    
    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        REPORT_STATUS();

        stats.update(align);

        if (!align.mapped || align.id != Standard::chrT)
        {
            return;
        }

        const Interval *match = nullptr;
        
        if (!align.spliced)
        {
            // Calculating statistics at the exon level
            match = matchExon(align, stats, fps);
        }
        else
        {
            // Calculating statistis at the intron level
            match = matchIntron(align, stats);
        }

        if (!match)
        {
            stats.unknowns.push_back(UnknownAlignment(align.id, align.l));
        }
    });

    o.info("Counting references");
    
    stats.pe.inferRefFromHist();
    stats.pi.inferRefFromHist();

    o.info("Calculating metrics for all sequins");

    /*
     * Calculating metrics for all sequins.
     */

    for (const auto &i : stats.eInters.data())
    {
        auto &m  = stats.sb.at(i.second.gID);
        auto &in = i.second;

        const auto &gID = i.second.gID;
        
        // Update the overall performance
        stats.pb.m.fp() += fps.at(gID);

        in.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                m.tp() += j - i;
                
                // Update the overall performance
                stats.pb.m.tp() += j - i;

                // Update the distribution
                stats.pb.h.at(gID)++;
            }
        });
        
        m.nr += in.l().length();
        m.nq += m.tp() + fps.at(i.second.gID);

        assert(m.nr >= m.tp());
        
        stats.pb.m.nr += in.l().length();
        stats.pb.m.nq  = stats.pb.m.tp() + stats.pb.m.fp();
    }

    o.info("Merging overlapping bases");

    assert(stats.pe.m.nr && stats.pi.m.nr && stats.pb.m.nr);

    /*
     * Calculate for the LOS
     */

    o.info("Calculating detection limit");

    stats.pe.hl = r.limitGene(stats.pe.h);
    stats.pi.hl = r.limitGene(stats.pi.h);
    stats.pb.hl = r.limitGene(stats.pb.h);

    stats.pe.m.sn();
    stats.pb.m.sn();
    stats.pi.m.sn();

    return stats;
}

void TAlign::report(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto stats = TAlign::stats(file, o);
    
    o.logInfo((boost::format("Exon: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pe.m.nr
                                          % stats.pe.m.nq
                                          % stats.pe.m.tp()
                                          % stats.pe.m.fp()
                                          % stats.pe.m.fn()
                                          % stats.pe.m.sn()
                                          % stats.pe.m.accuracy()).str());

    o.logInfo((boost::format("Intron: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pi.m.nr
                                          % stats.pi.m.nq
                                          % stats.pi.m.tp()
                                          % stats.pi.m.fp()
                                          % stats.pi.m.fn()
                                          % stats.pi.m.sn()
                                          % stats.pi.m.accuracy()).str());

    o.logInfo((boost::format("Base: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pb.m.nr
                                          % stats.pb.m.nq
                                          % stats.pb.m.tp()
                                          % stats.pb.m.fp()
                                          % stats.pb.m.fn()
                                          % stats.pb.m.sn()
                                          % stats.pb.m.accuracy()).str());

    /*
     * Write out summary statistics
     */

    {
        const auto summary = "Summary for dataset: %1%\n\n"
                             "   Unmapped:   %2% reads\n"
                             "   Experiment: %3% reads\n"
                             "   Synthetic:  %4% reads\n\n"
                             "   Reference:  %5% exons\n"
                             "   Reference:  %6% introns\n"
                             "   Reference:  %7% bases\n\n"
                             "   Query:      %8% exons\n"
                             "   Query:      %9% introns\n"
                             "   Query:      %10% bases\n\n"
                             "   Dilution:   %23%\n\n"
                             "   ***\n"
                             "   *** The following statistics are computed at the exon, intron and base level.\n"
                             "   ***\n"
                             "   ***   Exon level is defined as the performance per exon. An alignment that\n"
                             "   ***   is not mapped entirely withing an exon is considered as a FP. The\n"
                             "   ***   intron level is similar.\n"
                             "   ***\n"
                             "   ***   Base level is defined as the performance per nucleotide. A partial\n"
                             "   ***   mapped read will have FP and TP.\n"
                             "   ***\n\n"
                             "   -------------------- Exon level --------------------\n\n"
                             "   Sensitivity: %11%\n"
                             "   Accuarcy:    %12%\n"
                             "   Detection:   %13% (%14%)\n\n"
                             "   -------------------- Intron level --------------------\n\n"
                             "   Sensitivity: %15%\n"
                             "   Accuarcy:    %16%\n"
                             "   Detection:   %17% (%18%)\n\n"
                             "   -------------------- Base level --------------------\n\n"
                             "   Sensitivity: %19%\n"
                             "   Accuarcy:    %20%\n"
                             "   Detection:   %21% (%22%)\n";
        
        o.writer->open("TransAlign_summary.stats");
        o.writer->write((boost::format(summary) % file
                                                % stats.unmapped
                                                % stats.n_expT
                                                % stats.n_chrT
                                                % r.countSortedExons()
                                                % r.countSortedIntrons()
                                                % r.exonBase()
                                                % stats.pe.m.nq
                                                % stats.pi.m.nq
                                                % stats.pb.m.nq
                                                % stats.pe.m.sn()
                                                % stats.pe.m.accuracy()
                                                % stats.pe.hl.abund
                                                % stats.pe.hl.id
                                                % stats.pi.m.sn()
                                                % stats.pi.m.accuracy()
                                                % stats.pi.hl.abund
                                                % stats.pi.hl.id
                                                % stats.pb.m.sn()
                                                % stats.pb.m.accuracy()
                                                % stats.pb.hl.abund
                                                % stats.pb.hl.id
                                                % stats.dilution()).str());
        o.writer->close();
    }

    /*
     * Generating detailed statistics for each sequin
     */
    
    {
        o.writer->open("TransAlign_quins.stats");
        o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
        
        const auto format = "%1%\t%2%\t%3%\t%4%";

        o.writer->write((boost::format(format) % "ID"
                                               % "Exon"
                                               % "Intron"
                                               % "Base").str());
        
        for (const auto &i : stats.pe.h)
        {
            o.writer->write((boost::format(format) % i.first
                                                   % stats.pe.h.at(i.first)
                                                   % stats.pi.h.at(i.first)
                                                   % stats.pb.h.at(i.first)).str());
        }
        
        o.writer->close();
    }

    /*
     * Generating detailed logs for the histogram
     */
    
    {
        o.info("Generating detailed logs for the histogram");

        o.logInfo("\n\n--------------------- Exon Histogram ---------------------");
        for (const auto &i : stats.pe.h) { o.logInfo((boost::format("%1% %2%") % i.first % i.second).str()); }

        o.logInfo("\n\n--------------------- Intron Histogram ---------------------");
        for (const auto &i : stats.pi.h) { o.logInfo((boost::format("%1% %2%") % i.first % i.second).str()); }

        o.logInfo("\n\n--------------------- Base Histogram ---------------------");
        for (const auto &i : stats.pb.h) { o.logInfo((boost::format("%1% %2%") % i.first % i.second).str()); }
    }
}