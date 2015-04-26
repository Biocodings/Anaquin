#include "r_assembly.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "parsers/parser_gtf.hpp"

using namespace Spike;

template <typename Iter, typename F> static void extractIntrons(const Iter &exons, F f)
{
    Feature intr;
    
    for (auto i = 0; i < exons.size(); i++)
    {
        if (i)
        {
            intr = exons[i];
            intr.l = Locus(exons[i - 1].l.end, exons[i].l.start);
            
            // Intron is simply a non-transcribed region spliced between exons
            f(exons[i - 1], exons[i], intr);
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

    std::vector<Feature> q_exons;

    ParserGTF::parse(file, [&](const Feature &f)
    {
        switch (f.type)
        {
            case Exon:
            {
                // We'll need it to construct introns later
                q_exons.push_back(f);

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

                /*
                 * Classify at the base level
                 */
                
                stats.mb.nq() += f.l.length();
                stats.mb.tp() += std::min(f.l.length(), countOverlaps(s.r_exons, f));                
                stats.mb.fp()  = stats.mb.nq() - stats.mb.tp();

                break;
            }

            case Transcript:
            {
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

                /*
                 * Classify at the base level
                 */

                stats.mb.nq() += f.l.length();
                stats.mb.tp() += std::min(f.l.length(), countOverlaps_map(seqs, f));
                stats.mb.fp()  = stats.mb.nq() - stats.mb.tp();
                
                break;
            }

            default:
            {
                break;
            }
        }
    });

    assert(!s.r_introns.empty());

    // There is no guarantee that the exons in the query are sorted
    std::sort(q_exons.begin(), q_exons.end(), [&](const Feature &f1, const Feature &f2)
    {
        return f1.l.start < f2.l.start;
    });

    extractIntrons(q_exons, [&](const Feature &, const Feature &, Feature &i)
                   {
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
                       
                       /*
                        * Classify at the base level
                        */

                       stats.mb.nq() += i.l.length();
                       stats.mb.tp() += std::min(i.l.length(), countOverlaps(s.r_introns, i));
                       stats.mb.fp()  = stats.mb.nq() - stats.mb.tp();
                   });

    stats.mt.nr() = seqs.size();
    stats.me.nr() = s.r_exons.size();
    stats.mi.nr() = s.r_introns.size();
    stats.mb.nr() = 1000; // TODO: Need super-loci coding

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