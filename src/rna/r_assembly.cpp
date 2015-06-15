#include "rna/r_assembly.hpp"
#include "stats/expression.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Spike;

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

    // The structure depends on the mixture
    const auto seqs = s.r_sequin(options.mix);

    std::map<SequinID, std::vector<Feature>> q_exons_;
    std::vector<Feature> q_exons;

    Counts i = 0;

    ParserGTF::parse(file, [&](const Feature &f, const ParserProgress &)
    {
        // Don't bother unless the transcript is a sequin or it's been filtered
        if (!stats.ce.count(f.tID) || options.filters.count(f.tID))
        {
            return;
        }
        
        i++;

        switch (f.type)
        {
            case Exon:
            {
                q_exons.push_back(f);
                q_exons_[f.tID].push_back(f);

                /*
                 * Classify at the exon level
                 */

                if (classify(stats.pe.m, f, [&](const Feature &)
                {
                    return find(s.r_exons, f, Exact);
                }))
                {
                    stats.e_lc.at(f.l)++;
                    stats.ce.at(f.tID)++;
                }

                break;
            }

            case Transcript:
            {
                /*
                 * Classify at the transctipt level
                 */

                if (classify(stats.pt.m, f, [&](const Feature &)
                {
                    return find_map(seqs, f, Exact);
                }))
                {
                    stats.t_lc.at(f.tID)++;
                    stats.ct.at(f.tID)++;
                }

                break;
            }

            // There're many other possibilties in a GTF file, but we don't need those
            default: { break; }
        }
    });

    /*
     * Sort the exons in the query since there is no guarantee that those are sorted
     */
    
    for (auto &i : q_exons_)
    {
        CHECK_AND_SORT(i.second);
    }

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
                           stats.i_lc.at(i.l)++;
                           stats.ci[i.tID]++;
                       }
                   });
    
    /*
     * Classify at the base level
     */

    countBase(s.r_l_exons, q_exons, stats.pb.m, stats.cb);

    /*
     * Setting the known references
     */
    
    sums(stats.e_lc, stats.pe.m.nr);
    sums(stats.i_lc, stats.pi.m.nr);

    // The number of sequins is also the number of known transcripts
    stats.pt.m.nr = seqs.size();
    
    // The known base length is the total length of all known exons
    stats.pb.m.nr = s.r_c_exons;

    /*
     * Calculate for the LOS
     */

    stats.pe.s = Expression::analyze(stats.ce, seqs);
    stats.pt.s = Expression::analyze(stats.ct, seqs);
    stats.pb.s = Expression::analyze(stats.cb, s.r_pair(options.mix));
    stats.pi.s = Expression::analyze(stats.ci, seqs);

    /*
     * Write out statistics for various levels
     */

    AnalyzeReporter::report("assembly.base.stats", stats.pb, stats.cb, options.writer);
    AnalyzeReporter::report("assembly.exons.stats", stats.pe, stats.ce, options.writer);
    AnalyzeReporter::report("assembly.intron.stats", stats.pi, stats.ci, options.writer);
    AnalyzeReporter::report("assembly.transcripts.stats", stats.pt, stats.ct, options.writer);

    return stats;
}