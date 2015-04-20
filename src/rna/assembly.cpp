#include "assembly.hpp"
#include "classify.hpp"
#include "standard.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "parsers/parser_gtf.hpp"

#include <iostream>
using namespace Spike;

template <typename Iter, typename F> void extractIntrons(const Iter &exons, F f)
{
    Feature intr;
    
    for (auto i = 0; i < exons.size(); i++)
    {
        if (i)
        {
            intr = exons[i];
            intr.l = Locus(exons[i - 1].l.end, exons[i].l.start);
            f(exons[i - 1], exons[i], intr);
        }
    }
}

AssemblyStats Assembly::analyze(const std::string &file, const Assembly::Options &options)
{
    AssemblyStats stats;
    const auto &r = Standard::instance();

    auto cb = countsForSequins();
    auto ce = countsForSequins();
    auto ct = countsForSequins();
    auto ci = countsForSequins();
    
    std::vector<Feature> exons;
    
	ParserGTF::parse(file, [&](const Feature &f)
	{
        classify(stats, f, [&](const Feature &)
                 {
                     switch (f.type)
                     {
                         case Transcript:
                         {
                             assert(r.r_seqs_iA.count(f.iID));
                             const auto &seq = r.r_seqs_iA.at(f.iID);
                             
                             assert(cb.count(f.iID) && ct.count(f.iID));
                             
                             if (tfp(seq.l == f.l, &stats.mt))
                             {
                                 cb[f.iID]++;
                                 ct[f.iID]++;
                                 
                                 return true;
                             }
                             else
                             {
                                 return false;
                             }
                             
                             break;
                         }
                             
                         case Exon:
                         {
                             exons.push_back(f);
                             assert(cb.count(f.iID) && ce.count(f.iID));
                             
                             if (tfp(find(r.exons, f), &stats.me))
                             {
                                 cb[f.iID]++;
                                 ce[f.iID]++;
                                 
                                 return true;
                             }
                             else
                             {
                                 return false;
                             }
                             
                             break;
                         }
                             
                         default:
                         {
                             throw std::runtime_error("Unknown assembly type!");
                         }
                     }
                 });
    });

    assert(!r.introns.empty());
    
    extractIntrons(exons, [&](const Feature &, const Feature &, Feature &i)
    {
        classify(stats, i,
                 [&](const Feature &)
                 {
                     if (tfp(find(r.introns, i), &stats.mi))
                     {
                         assert(ci.count(i.iID));
                         ci[i.iID]++;
                         return true;
                     }
                     else
                     {
                         return false;
                     }
                 });
    });

    assert(stats.n && stats.nr + stats.nq == stats.n);

    
    const auto rb = Expression::analyze(cb);
    const auto re = Expression::analyze(ce);
    const auto rt = Expression::analyze(ct);
    const auto ri = Expression::analyze(ci);
    
    stats.sb = rb.sens(r.r_seqs_iA);
    stats.se = re.sens(r.r_seqs_iA);
    stats.st = rt.sens(r.r_seqs_iA);
    stats.si = ri.sens(r.r_seqs_iA);

    /*
     * Base-level statistics
     */
    
    AnalyzeReporter::reportSS("assembly.stats", stats, options.writer);
    
    /*
     * Exon-level statistics
     */

    /*
     * Transcript-level statistics
     */

    /*
     * Intron-level statistics
     */
    
    /*
     * Counting statistics
     */
    
    AnalyzeReporter::reportCounts("base.counts", cb, options.writer);
    AnalyzeReporter::reportCounts("exon.counts", ce, options.writer);
    AnalyzeReporter::reportCounts("transcripts.counts", ce, options.writer);
    AnalyzeReporter::reportCounts("introns.counts", ce, options.writer);

    return stats;
}