#include <fstream>
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

static std::string createFilteredGTF(const std::string &file)
{
    std::string line;
    const auto tmp = "ABCD.gtf"; //tmpnam(NULL);
    
    std::ofstream out(tmp);

    ParserGTF::parse(file, [&](const Feature &f, const std::string &l, const ParserProgress &)
    {
        if (f.id == Standard::instance().id)
        {
            out << l << std::endl;
        }
    });
    
    out.close();
    
    return tmp;
}

TAssembly::Stats TAssembly::report(const std::string &file, const Options &o)
{
    assert(!o.ref.empty() && !o.query.empty());

    /*
     * Comparing transcripts require constructing intron-chains, this is quite complicated.
     * We will reuse the code in cuffcompare. The idea is dirty but it works better than
     * than reinventing the wheel. However, we'll need to filter out only the features
     * belong to the synthetic chromosome.
     */

    o.logInfo("Creating a filtered transcript");

    const auto query = createFilteredGTF(file);
    o.logInfo("Filtered transcript: " + query + " has been created");

    o.logInfo("Invoking cuffcompare: " + o.ref);
    o.logInfo("Invoking cuffcompare: " + query);

    const int status = cuffcompare_main(o.ref.c_str(), o.query.c_str());

    //std::remove(query);
    
    if (status)
    {
        throw std::runtime_error("Failed to compare the given transcript. Please check the file and try again.");
    }
    
    TAssembly::Stats stats;
    const auto &r = Standard::instance().r_trans;

    std::vector<Feature> q_exons;
    std::map<SequinID, std::vector<Feature>> q_exons_;

    o.info("Parsing transcript");

    ParserGTF::parse(file, [&](const Feature &f, const std::string &, const ParserProgress &p)
    {
        if ((p.i % 1000000) == 0)
        {
            o.wait(std::to_string(p.i));
        }

        if (f.id != Standard::instance().id)
        {
            stats.n_hg38++;
            return;
        }
        else if ((p.i % 1000000) == 0)
        {
            o.wait(std::to_string(p.i));
        }
        
        stats.n_chrT++;

        switch (f.type)
        {
            /*
             * Classify at the exon level
             */

            case Exon:
            {
                const TransRef::ExonData *d;

                q_exons.push_back(f);
                q_exons_[f.tID].push_back(f);

                if (classify(stats.pe.m, f, [&](const Feature &)
                {
                    return (d = r.findExon(f.l, TransRef::Exact));
                }))
                {
                    stats.he.at(d->iID)++;
                }

                break;
            }

            /*
             * Classify at the transctipt level
             */

            case Transcript:
            {
                const TransData *match;

                if (classify(stats.pt.m, f, [&](const Feature &)
                {
                    return (match = r.seq(f.l));
                }))
                {
                    stats.ht.at(match->id)++;
                }

                break;
            }

            // There're many other possibilties in a GTF file, but we don't need them
            default: { break; }
        }
    });

    o.info("Generating introns");

    /*
     * Sort the query exons since there is no guarantee that those are sorted
     */

    for (auto &i : q_exons_)
    {
        CHECK_AND_SORT(i.second);
    }

    /*
     * Now that the query exons are sorted. We can extract and classify the introns for each pair
     * of successive exon.
     */
    
    const TransRef::IntronData *m;

    extractIntrons(q_exons_, [&](const Feature &, const Feature &, Feature &i)
    {
        if (classify(stats.pi.m, i, [&](const Feature &)
        {
            return (m  = r.findIntron(i.l, TransRef::Exact));
        }))
        {
            stats.hi.at(m->iID)++;
        }
    });

    o.info("Counting references");

    /*
     * Setting the known references
     */
    
    sums(stats.he, stats.pe.m.nr);
    sums(stats.hi, stats.pi.m.nr);

    /*
     * The counts for references for the transcript level is simply all the sequins.
     */

    stats.pt.m.nr = r.data().size();

    o.info("Merging overlapping bases");

    /*
     * The counts for query bases is the total non-overlapping length of all the exons in the experiment.
     * The number is expected to approach the reference length (calculated next) for a very large
     * experiment with sufficient coverage.
     */
    
    countBase(r.mergedExons(), q_exons, stats.pb.m, stats.hb);

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

    o.info("Calculating limit of sensitivity");
    
    stats.pe.s = r.limit(stats.he);
    stats.pt.s = r.limit(stats.ht);
    stats.pb.s = r.limitGene(stats.hb);
    stats.pi.s = r.limit(stats.hi);

    o.info("Generating statistics");

    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Genome: %2% reads\n"
                         "   Query: %3% reads\n"
                         "   Reference: %4% sequins\n\n"
                         "   Fuzzy: %5%\n\n"
                         "#--------------------|   Sn   |   Sp   |   Ss   |   fSn   |   fSp\n"
                         "    Exon level:       %6%     %7%     %8% (%9%)    %10%    %11%\n"
                         "    Intron level:       %12%     %13%     %14% (%15%)    %16%    %17%\n"
                         "    Base level:       %18%     %19%     %20% (%21%)    %22%    %23%\n"
                         "    Transcript level:       %24%     %25%     %26% (%27%)    %28%    %29%\n"
    ;
    
    o.writer->open("TransAssembly_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % "NA"
                                            % "NA"
                                            % r.data().size()
                                            % o.fuzzy
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
    o.writer->close();

    return stats;
}