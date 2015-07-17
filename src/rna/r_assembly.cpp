#include "rna/r_assembly.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Anaquin;

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

RAssembly::Stats RAssembly::analyze(const std::string &file, const Options &options)
{
    RAssembly::Stats stats;
    const auto &s = Standard::instance();

    // The sequins depends on the mixture
    const auto sequins = s.r_sequin(options.mix);

    std::vector<Feature> q_exons;
    std::map<SequinID, std::vector<Feature>> q_exons_;

    options.info("Parsing transcript");

    ParserGTF::parse(file, [&](const Feature &f, const ParserProgress &p)
    {
        if ((p.i % 1000000) == 0)
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
                const Feature *match;

                q_exons.push_back(f);
                q_exons_[f.tID].push_back(f);

                /*
                 * Classify at the exon level
                 */

                if (classify(stats.pe.m, f, [&](const Feature &)
                {
                    return (match = find(s.r_exons, f, Exact));
                }))
                {
                    stats.ce.at(match->tID)++;
                }
                else
                {
                    //options.logger->write((boost::format("[Exon]: %1% %2%") % std::to_string(f.l.start)
                      //                                                      % std::to_string(f.l.end)).str()) ;
                }

                break;
            }

            case Transcript:
            {
                const Sequin *match;

                /*
                 * Classify at the transctipt level
                 */

                if (classify(stats.pt.m, f, [&](const Feature &)
                {
                    return (match = find(sequins, f, Exact));
                }))
                {
                    stats.ct.at(match->id)++;
                }
                else
                {
                    options.logger->write((boost::format("[Transcript]: %1% %2%") % std::to_string(f.l.start)
                                                                                  % std::to_string(f.l.end)).str()) ;
                }

                break;
            }

            // There're many other possibilties in a GTF file, but we don't need those
            default: { break; }
        }
    });

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
                         return find(s.r_introns, i, Exact);
                     }))
        {
            stats.ci[i.tID]++;
        }
    });

    options.info("Counting references");

    /*
     * Setting the known references
     */
    
    sums(stats.ce, stats.pe.m.nr);
    sums(stats.ci, stats.pi.m.nr);

    /*
     * The counts for references for the transcript level is simply all the sequins.
     */

    stats.pt.m.nr = sequins.size();

    options.info("Merging overlapping bases");

    /*
     * The counts for query bases is the total non-overlapping length of all the exons in the experiment.
     * The number is expected to approach the reference length (calculated next) for a very large
     * experiment with sufficient coverage.
     */
    
    countBase(s.r_l_exons, q_exons, stats.pb.m, stats.cb);

    /*
     * The counts for references is the total length of all known non-overlapping exons.
     * For example, if we have the following exons:
     *
     *    {1,10}, {50,55}, {70,74}
     *
     * The length of all the bases is 10+5+4 = 19.
     */

    stats.pb.m.nr = s.r_c_exons;

    /*
     * Calculate for the LOS
     */

    options.info("Calculating LOS - exon level");
    stats.pe.s = Expression::analyze(stats.ce, sequins);
    
    options.info("Calculating LOS - transcript level");
    stats.pt.s = Expression::analyze(stats.ct, sequins);
    
    options.info("Calculating LOS - base level");
    stats.pb.s = Expression::analyze(stats.cb, s.r_gene(options.mix));
    
    options.info("Calculating LOS - intron level");
    stats.pi.s = Expression::analyze(stats.ci, sequins);

    options.info("Generating base statistics");
    AnalyzeReporter::stats("rna_assembly.base.stats", stats.pb, stats.cb, options.writer);

    options.info("Generating exon statistics");
    AnalyzeReporter::stats("rna_assembly.exons.stats", stats.pe, stats.ce, options.writer);
    
    options.info("Generating intron statistics");
    AnalyzeReporter::stats("rna_assembly.intron.stats", stats.pi, stats.ci, options.writer);
    
    options.info("Generating transcript statistics");
    AnalyzeReporter::stats("rna_assembly.transcripts.stats", stats.pt, stats.ct, options.writer);

    return stats;
}