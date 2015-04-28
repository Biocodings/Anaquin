#include "r_assembly.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "parsers/parser_gtf.hpp"
#include <iostream>
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
                ir = ts.second[i];
                ir.l = Locus(ts.second[i - 1].l.end, ts.second[i].l.start);

                if (ts.second[i-1].iID == ts.second[i].iID)
                {
                    // Intron is simply a non-transcribed region spliced between exons
                    f(ts.second[i-1], ts.second[i], ir);
                }
            }
        }
    }
}

RAssemblyStats RAssembly::analyze(const std::string &file, const Options &options)
{
    RAssemblyStats stats;
    const auto &s = Standard::instance();

    auto cb = RAnalyzer::counter(Isoform, options.mix);
    auto ce = RAnalyzer::counter(Isoform, options.mix);
    auto ct = RAnalyzer::counter(Isoform, options.mix);
    auto ci = RAnalyzer::counter(Isoform, options.mix);

    // The structure depends on the mixture
    const auto seqs = s.r_sequin(options.mix);

    std::map<SequinID, std::vector<Feature>> q_exons_;
    std::vector<Feature> q_exons, q_trans, q_introns;

    ParserGTF::parse(file, [&](const Feature &f)
    {
        if (options.filters.count(f.iID))
        {
            return;
        }

        switch (f.type)
        {
            case Exon:
            {
                q_exons.push_back(f);
                q_exons_[f.iID].push_back(f);

                /*
                 * Classify at the exon level
                 */
                
                if (classify(stats.me, f, [&](const Feature &)
                {
                    return find(s.r_exons, f, ExactRule);
                }))
                {
                    ce.at(f.iID)++;
                }

                assert(stats.mb.nq() >= stats.mb.tp());
                break;
            }

            case Transcript:
            {
                q_trans.push_back(f);
                
                /*
                 * Classify at the transctipt level
                 */

                if (classify(stats.mt, f, [&](const Feature &)
                {
                    return find_map(seqs, f, ExactRule);
                }))
                {
                    ct.at(f.iID)++;
                }

                break;
            }

            default:
            {
                break;
            }
        }
    });

    assert(!s.r_introns.empty());

    /*
     * Sort the exons in the query since there is no guarantee that those are sorted.
     */
    
    for (auto &i : q_exons_)
    {
        std::sort(i.second.begin(), i.second.end(), [&](const Feature &f1, const Feature &f2)
        {
            return f1.l.start < f2.l.start;
        });
    }

    extractIntrons(q_exons_, [&](const Feature &, const Feature &, Feature &i)
                   {
                       q_introns.push_back(i);

                       /*
                        * Classify at the intron level
                        */
                       
                       if (classify(stats.mi, i, [&](const Feature &)
                       {
                           return find(s.r_introns, i, ExactRule);
                       }))
                       {
                           ci[i.iID]++;
                       }
                   });
    
    /*
     * Classify at the base level
     */

    countBase(s.r_l_exons,   q_exons,   stats.mb);
    countBase(s.r_l_trans,   q_trans,   stats.mb);
    countBase(s.r_l_introns, q_introns, stats.mb);

    /*
     * Setting the known references
     */

    stats.mt.nr() = seqs.size();
    stats.me.nr() = s.r_exons.size();
    stats.mi.nr() = s.r_introns.size();
    stats.mb.nr() = s.r_c_exons + s.r_c_trans + s.r_c_introns;

    assert(stats.mb.nr() >= stats.mb.tp());
    assert(stats.mt.nr() >= stats.mt.tp());
    assert(stats.me.nr() >= stats.me.tp());
    assert(stats.mi.nr() >= stats.mi.tp());

    assert(stats.me.nq() == q_exons.size());
    assert(stats.me.nq() >= stats.me.tp());
    assert(s.r_exons.size() >= stats.me.tp());

    assert(stats.mt.nq() >= stats.mt.tp());
    assert(stats.me.nq() >= stats.me.tp());
    assert(stats.mi.nq() >= stats.mi.tp());

    /*
     * Calculate for the sensitivity
     */

    stats.se = Expression::analyze(ce, seqs);
    stats.st = Expression::analyze(ct, seqs);
    stats.sb = Expression::analyze(cb, seqs);
    stats.si = Expression::analyze(ci, seqs);

    /*
     * Report for the statistics
     */
    
    const auto &writer = options.writer;

    // Report the for base level
    AnalyzeReporter::reportClassify("assembly.base.stats", stats.mb, stats.sb, cb, writer);

    // Report the for exons level
    AnalyzeReporter::reportClassify("assembly.exons.stats", stats.me, stats.se, ce, writer);

    // Report the for transcripts level
    AnalyzeReporter::reportClassify("assembly.transcripts.stats", stats.mt, stats.st, ct, writer);

    // Report the for intron level
    AnalyzeReporter::reportClassify("assembly.intron.stats", stats.mi, stats.si, ci, writer);

    return stats;
}