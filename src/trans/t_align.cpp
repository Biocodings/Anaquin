#include "trans/t_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

TAlign::Stats TAlign::report(const FileName &file, const Options &o)
{
    TAlign::Stats stats;
    const auto &r = Standard::instance().r_trans;

    std::vector<Alignment> exons, introns;

    o.analyze(file);

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (!align.i)
        {
            if      (!align.mapped)                       { stats.unmapped++; }
            else if (align.id != Standard::instance().id) { stats.n_expT++;   }
            else                                          { stats.n_chrT++;   }
        }

        if (!align.mapped || align.id != Standard::instance().id)
        {
            return;
        }
        
        /*
         * Collect statistics at the exon level
         */

        if (!align.spliced)
        {
            exons.push_back(align);

            const TransRef::ExonData *d;

            if (classify(stats.pe.m, align, [&](const Alignment &)
            {
                return (d = r.findExon(align.l, Contains));
            }))
            {
                stats.he.at(d->gID)++;
            }
        }

        /*
         * Collect statistics at the intron level
         */

        else
        {
            introns.push_back(align);

            const TransRef::IntronData *d;

            if (classify(stats.pi.m, align, [&](const Alignment &)
            {
                return (d = r.findIntron(align.l, Exact));
            }))
            {
                stats.hi.at(d->gID)++;
            }
        }
    });

    o.info("Counting references");

    /*
     * Calculate for references. The idea is similar to cuffcompare, each true-positive is counted
     * as a reference. Anything that is undetected in the experiment will be counted as a single reference.
     */

    sums(stats.he, stats.pe.m.nr);
    sums(stats.hi, stats.pi.m.nr);
    
    o.info("Merging overlapping bases");

    /*
     * Counts at the base-level is the non-overlapping region of all the exons
     */

    countBase(r.mergedExons(), exons, stats.pb.m, stats.hb);

    /*
     * The counts for references is the total length of all known non-overlapping exons.
     * For example, if we have the following exons:
     *
     *    {1,10}, {50,55}, {70,74}
     *
     * The length of all the bases is 10+5+4 = 19.
     */

    stats.pb.m.nr = r.exonBase();
    assert(stats.pe.m.nr && stats.pi.m.nr && stats.pb.m.nr);

    /*
     * Calculate for the LOS
     */

    o.info("Calculating detection limit");

    stats.pe.s = r.limitGene(stats.he);
    stats.pi.s = r.limitGene(stats.hi);
    stats.pb.s = r.limitGene(stats.hb);

    assert(stats.he.size() == stats.hi.size());
    assert(stats.hi.size() == stats.hb.size());

    stats.pe.m.sn();
    stats.pb.m.sn();
    stats.pi.m.sn();
    
    o.logInfo((boost::format("Exon: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pe.m.nr
                                          % stats.pe.m.nq
                                          % stats.pe.m.tp()
                                          % stats.pe.m.fp()
                                          % stats.pe.m.fn()
                                          % stats.pe.m.sn()
                                          % stats.pe.m.sp()).str());

    o.logInfo((boost::format("Base: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pb.m.nr
                                          % stats.pb.m.nq
                                          % stats.pb.m.tp()
                                          % stats.pb.m.fp()
                                          % stats.pb.m.fn()
                                          % stats.pb.m.sn()
                                          % stats.pb.m.sp()).str());
    
    o.logInfo((boost::format("Intron: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pi.m.nr
                                          % stats.pi.m.nq
                                          % stats.pi.m.tp()
                                          % stats.pi.m.fp()
                                          % stats.pi.m.fn()
                                          % stats.pi.m.sn()
                                          % stats.pi.m.sp()).str());

    /*
     * Write out summary statistics
     */

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Unmapped:  %2% reads\n"
                         "   Genome:    %3% reads\n"
                         "   Synthetic: %4% reads\n\n"
                         "   Reference: %5% exons\n"
                         "   Reference: %6% introns\n"
                         "   Reference: %7% bases\n\n"
                         "   Query: %8% exons\n"
                         "   Query: %9% introns\n"
                         "   Query: %10% bases\n\n"
                         "#--------------------|   Sn   |  Sp  |  Ss  \n"
                         "    Exon level:\t%11%\t%12%\t%13% (%14%)\n"
                         "    Intron level:\t%15%\t%16%\t%17% (%18%)\n"
                         "    Base level:\t%19%\t%20%\t%21% (%22%)\n"
                         "\n"
                         "Dilution:\t\t%23%\n"
    ;

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
                                            % stats.pe.m.sp()
                                            % stats.pe.s.abund
                                            % stats.pe.s.id
                                            % stats.pi.m.sn()
                                            % stats.pi.m.sp()
                                            % stats.pi.s.abund
                                            % stats.pi.s.id
                                            % stats.pb.m.sn()
                                            % stats.pb.m.sp()
                                            % stats.pb.s.abund
                                            % stats.pb.s.id
                                            % stats.dilution()
                                        ).str());
    o.writer->close();

    /*
     * Generating detailed statistics for each sequin
     */
    
    o.writer->open("TransAlign_quins.stats");
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());

    const auto format = "%1%\t%2%\t%3%\t%4%";
    o.writer->write((boost::format(format) % "ID" % "Exon" % "Intron" % "Base").str());

    for (const auto &i : stats.he)
    {
        o.writer->write((boost::format(format) % i.first
                                               % stats.he.at(i.first)
                                               % stats.hi.at(i.first)
                                               % stats.hb.at(i.first)).str());
    }
    
    o.writer->close();

    /*
     * Generating detailed logs for the histogram
     */
    
    {
        o.info("Generating detailed logs for the histogram");

        o.logInfo("\n\n--------------------- Exon Histogram ---------------------");
        for (const auto &i : stats.he) { o.logInfo((boost::format("%1% %2%") % i.first % i.second).str()); }

        o.logInfo("\n\n--------------------- Intron Histogram ---------------------");
        for (const auto &i : stats.hi) { o.logInfo((boost::format("%1% %2%") % i.first % i.second).str()); }

        o.logInfo("\n\n--------------------- Base Histogram ---------------------");
        for (const auto &i : stats.hb) { o.logInfo((boost::format("%1% %2%") % i.first % i.second).str()); }
    }

	return stats;
}