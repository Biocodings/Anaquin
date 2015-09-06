#include "trans/t_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

// Find the matching intron by locus given a spliced alignment
static bool findIntron(const Alignment &align, Anaquin::TransRef::IntronData &d)
{
    assert(align.spliced);
    const auto &s = Standard::instance();

    for (const auto &i : s.r_trans.sortedIntrons())
    {
        if (align.l == i.l)
        {
            d = i;
            return true;
        }
    }
    
    return false;
}

TAlign::Stats TAlign::analyze(const std::string &file, const Options &o)
{
    TAlign::Stats stats;
    const auto &s = Standard::instance();

    std::vector<Alignment> exons, introns;

    o.info("Parsing alignment file");

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if (!align.i && (p.i % 1000000) == 0)
        {
            o.wait(std::to_string(p.i));
        }
        
        if (align.id != s.id && !align.i)
        {
            stats.n_genome++;
        }

        if (!align.mapped || align.id != s.id)
        {
            return;
        }
        else if (!align.i)
        {
            stats.n_chrT++;            
        }
        
        // Whether the read has mapped to a feature correctly
        bool succeed = false;

        o.logInfo((boost::format("%1% %2% %3%") % align.id % align.l.start % align.l.end).str());

        /*
         * Collect statistics at the exon level
         */

        if (!align.spliced)
        {
            exons.push_back(align);

            const TransRef::ExonData *d;

            if (classify(stats.pe.m, align, [&](const Alignment &)
            {
                return (d = s.r_trans.findExon(align.l));
            }))
            {
                stats.he.at(s.seq2base.at(d->iID))++;
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
                return (d = s.r_trans.findIntron(align.l));
            }))
            {
                stats.hi.at(s.seq2base.at(d->iID))++;
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

    // TODOcountBase(s.r_l_exons, exons, stats.pb.m, stats.hb);

    /*
     * The counts for references is the total length of all known non-overlapping exons.
     * For example, if we have the following exons:
     *
     *    {1,10}, {50,55}, {70,74}
     *
     * The length of all the bases is 10+5+4 = 19.
     */
    
    // TODOstats.pb.m.nr = s.r_c_exons;

    assert(stats.pe.m.nr && stats.pi.m.nr && stats.pb.m.nr);

    /*
     * Calculate for the LOS
     */

    o.info("Calculating limit of sensitivity");

    stats.pe.s = Expression::analyze(stats.he, s.bases_1);
    stats.pi.s = Expression::analyze(stats.hi, s.bases_1);
    stats.pb.s = Expression::analyze(stats.hb, s.bases_1);

    /*
     * Write out summary statistics
     */

    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Genome: %2% reads\n"
                         "   Query: %3% reads\n"
                         "Reference: %4% exons\n\n"
                         "Fuzzy: %5%\n\n"
                         "#--------------------|   Sn   |  Sp   |  fSn |  fSp\n"
                         "    Exon level:       %6%     %7%     %8%    %9%\n"
                         "    Intron level:       %10%     %11%     %12%    %13%\n"
                         "    Base level:       %14%     %15%     %16%    %17%\n"
                         "\n"
                         "Dilution:     %18%\n"
                         "Detection Limit:     %19% (%20%)\n"
    ;

    o.writer->open("TransAlign_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.n_genome
                                            % stats.n_chrT
                                            % stats.pe.m.nr
                                            % o.fuzzy
                                            % stats.pe.m.sn()
                                            % stats.pe.m.sp()
                                            % "NA"
                                            % "NA"
                                            % stats.pi.m.sn()
                                            % stats.pi.m.sp()
                                            % "NA"
                                            % "NA"
                                            % stats.pb.m.sn()
                                            % stats.pb.m.sp()
                                            % "NA"
                                            % "NA"
                                            % stats.dilution()
                                            % stats.pe.s.abund
                                            % stats.pe.s.id
                                        ).str());
    o.writer->close();

	return stats;
}