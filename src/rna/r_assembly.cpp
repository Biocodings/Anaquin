#include "classify.hpp"
#include "standard.hpp"
#include "r_assembly.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "parsers/parser_gtf.hpp"

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
            
            // Intron is simply a non-transcribed region spliced between exons
            f(exons[i - 1], exons[i], intr);
        }
    }
}

RAssemblyStats RAssembly::analyze(const std::string &file, const Options &options)
{
    RAssemblyStats stats;
    const auto &r = Standard::instance();

    auto cb = countsForSequins();
    auto ce = countsForSequins();
    auto ct = countsForSequins();
    auto ci = countsForSequins();

    // The structure depends on the mixture
    const auto seq = r.r_mix_sequin(options.mix);

    std::vector<Feature> exons;

    

    ParserGTF::parse(file, [&](const Feature &f)
    {

    });
    
    
    
    
    
    
	ParserGTF::parse(file, [&](const Feature &f)
	{
        classify(stats, f, [&](const Feature &)
                 {
                     switch (f.type)
                     {
                         case Transcript:
                         {
                             assert(seq.count(f.iID));
                             const auto &s = seq.at(f.iID);
                             
                             assert(cb.count(f.iID) && ct.count(f.iID));

                             if (tfp(s.l == f.l, &stats.mt))
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

    assert(stats.n() && stats.nr + stats.nq == stats.n());

    const auto rb = Expression::analyze(cb);
    const auto re = Expression::analyze(ce);
    const auto rt = Expression::analyze(ct);
    const auto ri = Expression::analyze(ci);
    
    stats.sb = rb.sens(r.r_seqs_iA);
    stats.se = re.sens(r.r_seqs_iA);
    stats.st = rt.sens(r.r_seqs_iA);
    stats.si = ri.sens(r.r_seqs_iA);

    const auto &writer = options.writer;

    // Report the for base level
    AnalyzeReporter::reportClassify("assembly.base.stats", stats.dilution(), stats.mb, stats.sb, cb, writer);

    // Report the for exons level
    AnalyzeReporter::reportClassify("assembly.exons.stats", stats.dilution(), stats.me, stats.se, ce, writer);

    // Report the for transcripts level
    AnalyzeReporter::reportClassify("assembly.transcripts.stats", stats.dilution(), stats.mt, stats.st, ct, writer);

    // Report the for intron level
    AnalyzeReporter::reportClassify("assembly.intron.stats", stats.dilution(), stats.mi, stats.si, ci, writer);
    
    return stats;
}