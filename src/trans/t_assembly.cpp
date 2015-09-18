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
    const auto tmp = tmpnam(NULL);
    
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

    o.info("Generating a filtered transcript");

    const auto query = createFilteredGTF(file);
    o.logInfo("Filtered transcript: " + query + " has been created");

    o.logInfo("Invoking cuffcompare: " + o.ref);
    o.logInfo("Invoking cuffcompare: " + query);

    const int status = cuffcompare_main(o.ref.c_str(), query.c_str());

    if (status)
    {
        throw std::runtime_error("Failed to assess the given transcript. Please check the file and try again.");
    }
    
    TAssembly::Stats stats;
    const auto &r = Standard::instance().r_trans;

    std::vector<Feature> q_exons;
    std::map<SequinID, std::vector<Feature>> q_exons_;

    o.info("Parsing transcript");

    Confusion t;
    
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
        
        stats.n_chrT++;

        switch (f.type)
        {
            /*
             * Classify at the exon level
             */

            case Exon:
            {
                const TransRef::ExonData *match;

                q_exons.push_back(f);
                q_exons_[f.tID].push_back(f);

                if (classify(t, f, [&](const Feature &)
                {
                    return (match = r.findExon(f.l, TransRef::Exact));
                }))
                {
                    stats.he.at(match->iID)++;
                }

                break;
            }

            /*
             * Classify at the transctipt level
             */

            case Transcript:
            {
                const TransData *match;

                if (classify(t, f, [&](const Feature &)
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
    
    const TransRef::IntronData *match;

    extractIntrons(q_exons_, [&](const Feature &, const Feature &, Feature &i)
    {
        if (classify(t, i, [&](const Feature &)
        {
            return (match  = r.findIntron(i.l, TransRef::Exact));
        }))
        {
            stats.hi.at(match->iID)++;
        }
    });

    /*
     * Calculate for the LOS
     */

    o.info("Calculating limit of sensitivity");
    
    stats.se = r.limit(stats.he);
    stats.st = r.limit(stats.ht);
    stats.sb = r.limitGene(stats.hb);
    stats.si = r.limit(stats.hi);

    o.info("Generating statistics");

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Genome: %2% features\n"
                         "   Query: %3% features\n\n"
                         "   Reference: %4% exons\n"
                         "   Reference: %5% introns\n\n"
                         "#--------------------|   Sn   |   Sp   |   Ss \n"
                         "    Exon level:\t%6%\t%7%\t%8% (%9%)\n"
                         "    Intron level:\t%10%\t%11%\t%12% (%13%)\n"
                         "    Base level:\t%14%\t%15%\t%16% (%17%)\n"
                         "    Transcript level:\t%18%\t%19%\t%20% (%21%)\n"
    ;

    o.writer->open("TransAssembly_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.n_hg38
                                            % stats.n_chrT
                                            % r.data().size()
                                            % r.countSortedIntrons()
                                            % (__cmp__.e_sp / 100.0)
                                            % (__cmp__.e_sn / 100.0)
                                            % (stats.se.id.empty() ? "-" : std::to_string(stats.se.abund))
                                            % stats.se.id
                                            % (__cmp__.i_sp / 100.0)
                                            % (__cmp__.i_sn / 100.0)
                                            % (stats.si.id.empty() ? "-" : std::to_string(stats.si.abund))
                                            % stats.si.id
                                            % (__cmp__.b_sp / 100.0)
                                            % (__cmp__.b_sn / 100.0)
                                            % (stats.sb.id.empty() ? "-" : std::to_string(stats.sb.abund))
                                            % stats.sb.id
                                            % (__cmp__.t_sp / 100.0)
                                            % (__cmp__.t_sn / 100.0)).str());
    o.writer->close();

    return stats;
}