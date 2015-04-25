#include "classify.hpp"
#include "standard.hpp"
#include "r_assembly.hpp"
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

    std::vector<Feature> exons;

    ParserGTF::parse(file, [&](const Feature &f)
    {
        switch (f.type)
        {
            case Exon:
            {
                // We'll need it to construct introns later
                exons.push_back(f);
                
                /*
                 * Classify at the exon level
                 */
                
                if (classify(stats.me, f, [&](const Feature &)
                {
                    return find(s.r_exons, f);
                }))
                {
                    ce.at(f.iID)++;
                }
                
                break;
            }

            case Transcript:
            {
                /*
                 * Classify at the transctipt level
                 */
                
                if (classify(stats.mt, f, [&](const Feature &)
                {
                    return find_map(seqs, f);
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

        /*
         * Classify at the base level. The
         */
        
        
    });

    stats.me.nr() = s.r_exons.size();
    stats.mt.nr() = seqs.size();

    assert(stats.me.nq() == exons.size());
    assert(stats.me.nq() >= stats.me.tp());
    assert(s.r_exons.size() >= stats.me.tp());

    assert(stats.mt.nq() >= stats.mt.tp());

    /*
     * Calculate for the sensitivity
     */

    stats.se = Expression::analyze(ce, seqs);
    stats.st = Expression::analyze(ct, seqs);

    /*
     * Report for the statistics
     */
    
    const auto &writer = options.writer;

    // Report the for base level
    //AnalyzeReporter::reportClassify("assembly.base.stats", stats.dilution(), stats.mb, stats.sb, cb, writer);
    
    // Report the for exons level
    AnalyzeReporter::reportClassify("assembly.exons.stats", stats.me, stats.se, ce, writer);

    // Report the for transcripts level
    AnalyzeReporter::reportClassify("assembly.transcripts.stats", stats.mt, stats.st, ct, writer);

    // Report the for intron level
    //AnalyzeReporter::reportClassify("assembly.intron.stats", stats.dilution(), stats.mi, stats.si, ci, writer);

    
    
    
//	ParserGTF::parse(file, [&](const Feature &f)
//	{
//        classify(stats, f, [&](const Feature &)
//                 {
//                     switch (f.type)
//                     {
//                         case Transcript:
//                         {
//                             assert(seq.count(f.iID));
//                             const auto &s = seq.at(f.iID);
//                             
//                             assert(cb.count(f.iID) && ct.count(f.iID));
//
//                             if (tfp(s.l == f.l, &stats.mt))
//                             {
//                                 cb[f.iID]++;
//                                 ct[f.iID]++;
//                                 
//                                 return true;
//                             }
//                             else
//                             {
//                                 return false;
//                             }
//
//                             break;
//                         }
//                             
//                         case Exon:
//                         {
//                             break;
//                         }
//                             
//                         default:
//                         {
//                             throw std::runtime_error("Unknown assembly type!");
//                         }
//                     }
//                 });
//    });

    
    
    
   // assert(!r.introns.empty());
    
    extractIntrons(exons, [&](const Feature &, const Feature &, Feature &i)
    {
//        classify(stats, i,
//                 [&](const Feature &)
//                 {
//                     if (tfp(find(r.introns, i), &stats.mi))
//                     {
//                         assert(ci.count(i.iID));
//                         ci[i.iID]++;
//                         return true;
//                     }
//                     else
//                     {
//                         return false;
//                     }
//                 });
    });

    //assert(stats.n() && stats.nr + stats.nq == stats.n());
/*
    const auto rb = Expression::analyze(cb);
    const auto ri = Expression::analyze(ci);

//    stats.sb = rb.sens(r.r_seqs_iA);
    //stats.si = ri.sens(r.r_seqs_iA);
*/
    return stats;
}