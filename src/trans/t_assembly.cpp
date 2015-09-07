#include "data/compare.hpp"
#include "trans/t_assembly.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Anaquin;

// Defined for cuffcompare
Compare __cmp__;

// Defined for cuffcompare
extern int cuffcompare_main(const char *ref, const char *query);

template <typename F> static void extractIntrons(const std::map<SequinID, std::vector<Feature>> &x, F f)
{
    Feature ir;

    for (const auto & ts : x)
    {
        for (auto i = 0; i < ts.second.size(); i++)
        {
            if (i)
            {
                if (ts.second[i-1].tID == ts.second[i].tID)
                {
                    ir = ts.second[i];
                    ir.l = Locus(ts.second[i - 1].l.end + 1, ts.second[i].l.start - 1);
                    f(ts.second[i-1], ts.second[i], ir);
                }
            }
        }
    }
}

TAssembly::Stats TAssembly::analyze(const std::string &file, const Options &options)
{
    assert(!options.ref.empty() && !options.query.empty());
    
    /*
     * Comparing transcripts require constructing intron-chains, this is quite complicated.
     * We will reuse the code in cuffcompare. The implementation is dirty but it works better
     * than reinventing the wheel.
     */
    
    options.logInfo("Invoking cuffcompare: " + options.ref);
    options.logInfo("Invoking cuffcompare: " + options.query);

    const int status = cuffcompare_main(options.ref.c_str(), options.query.c_str());

    if (status)
    {
        throw std::runtime_error("Failed to compare the given transcript. Please check the file and try again.");
    }
    
    TAssembly::Stats stats;
    const auto &r = Standard::instance().r_trans;

    std::vector<Feature> q_exons;
    std::map<SequinID, std::vector<Feature>> q_exons_;

    options.info("Parsing transcript");

    ParserGTF::parse(file, [&](const Feature &f, const ParserProgress &p)
    {
        if (f.id != Standard::instance().id)
        {
            return;
        }
        else if ((p.i % 1000000) == 0)
        {
            options.wait(std::to_string(p.i));
        }

        // Don't bother unless the transcript is a sequin or it's been filtered
        if (options.filters.count(f.tID))
        {
            return;
        }
        
        options.logInfo((boost::format("%1% %2% %3%") % p.i % f.l.start % f.l.end).str());

        switch (f.type)
        {
            case Exon:
            {
                const TransRef::ExonData *d;

                q_exons.push_back(f);
                q_exons_[f.tID].push_back(f);

                /*
                 * Classify at the exon level
                 */

                if (classify(stats.pe.m, f, [&](const Feature &)
                {
                    return (d = r.findExon(f.l, TransRef::Exact));
                }))
                {
                    stats.he.at(d->iID)++;
                }

                break;
            }

            case Transcript:
            {
                const TransData *match;

                /*
                 * Classify at the transctipt level
                 */

                if (classify(stats.pt.m, f, [&](const Feature &)
                {
                    //return (match = find(sequins, f, Exact));
                    return (match = r.seq(f.l));
                }))
                {
                    stats.ht.at(match->id)++;
                }
                else
                {
                    options.logger->write(
                        (boost::format("[Transcript]: %1% %2%") % std::to_string(f.l.start)
                                                                % std::to_string(f.l.end)).str()) ;
                }

                break;
            }

            // There're many other possibilties in a GTF file, but we don't need those
            default: { break; }
        }
    });

    options.info("Transcript parsed");

    /*
     * Sort the query exons since there is no guarantee that those are sorted
     */

    for (auto &i : q_exons_)
    {
        CHECK_AND_SORT(i.second);
    }

    /*
     * Now that the query exons are sorted. We can extract the introns between each pair
     * of successive exon.
     */

    extractIntrons(q_exons_, [&](const Feature &, const Feature &, Feature &i)
    {
        /*
         * Classify at the intron level
         */
        
        if (classify(stats.pi.m, i, [&](const Feature &)
        {
            return r.findIntron(i.l, TransRef::Exact);
        }))
        {
            stats.hi[i.tID]++;
        }
    });

    options.info("Counting references");

    /*
     * Setting the known references
     */
    
    sums(stats.he, stats.pe.m.nr);
    sums(stats.hi, stats.pi.m.nr);

    /*
     * The counts for references for the transcript level is simply all the sequins.
     */

    stats.pt.m.nr = r.data().size();

    options.info("Merging overlapping bases");

    /*
     * The counts for query bases is the total non-overlapping length of all the exons in the experiment.
     * The number is expected to approach the reference length (calculated next) for a very large
     * experiment with sufficient coverage.
     */
    
    //countBase(s.r_l_exons, q_exons, stats.pb.m, stats.hb);

    /*
     * The counts for references is the total length of all known non-overlapping exons.
     * For example, if we have the following exons:
     *
     *    {1,10}, {50,55}, {70,74}
     *
     * The length of all the bases is 10+5+4 = 19.
     */

    //stats.pb.m.nr = s.r_c_exons;

    /*
     * Calculate for the LOS
     */

    options.info("Calculating limit of sensitivity");
    
    stats.pe.s = Expression_::calculate(stats.he, r);
    stats.pt.s = Expression_::calculate(stats.ht, r);
    stats.pb.s = Expression_::calculate(stats.hb, r);
    stats.pi.s = Expression_::calculate(stats.hi, r);

    options.info("Generating statistics");

    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Genome: %2% reads\n"
                         "   Query: %3% reads\n"
                         "Reference: %4% sequins\n\n"
                         "Fuzzy: %5%\n\n"
                         "#--------------------|   Sn   |   Sp   |   Ss   |   fSn   |   fSp\n"
                         "    Exon level:       %6%     %7%     %8% (%9%)    %10%    %11%\n"
                         "    Intron level:       %12%     %13%     %14% (%15%)    %16%    %17%\n"
                         "    Base level:       %18%     %19%     %20% (%21%)    %22%    %23%\n"
                         "    Transcript level:       %24%     %25%     %26% (%27%)    %28%    %29%\n"
    ;
    
    options.writer->open("TransAssembly_summary.stats");
    options.writer->write((boost::format(summary) % file
                                                  % "NA"
                                                  % "NA"
                                                  % r.data().size()
                                                  % options.fuzzy
                                                  % (__cmp__.e_sp / 100.0)
                                                  % (__cmp__.e_sn / 100.0)
                                                  % (stats.pe.s.id.empty() ? "-" : std::to_string(stats.pe.s.abund))
                                                  % stats.pe.s.id
                                                  % (__cmp__.e_fsp / 100.0)
                                                  % (__cmp__.e_fsn / 100.0)
                                                  % (__cmp__.i_sp / 100.0)
                                                  % (__cmp__.i_sn / 100.0)
                                                  % (stats.pi.s.id.empty() ? "-" : std::to_string(stats.pi.s.abund))
                                                  % stats.pi.s.id
                                                  % (__cmp__.i_fsp / 100.0)
                                                  % (__cmp__.i_fsn / 100.0)
                                                  % (__cmp__.b_sp / 100.0)
                                                  % (__cmp__.b_sn / 100.0)
                                                  % (stats.pb.s.id.empty() ? "-" : std::to_string(stats.pb.s.abund))
                                                  % stats.pb.s.id
                                                  % "-"
                                                  % "-"
                                                  % (__cmp__.t_sp / 100.0)
                                                  % (__cmp__.t_sn / 100.0)
                                                  % (stats.pt.s.id.empty() ? "-" : std::to_string(stats.pt.s.abund))
                                                  % stats.pt.s.id
                                                  % (__cmp__.t_fsp / 100.0)
                                                  % (__cmp__.t_fsn / 100.0)).str());
    options.writer->close();

    return stats;
}