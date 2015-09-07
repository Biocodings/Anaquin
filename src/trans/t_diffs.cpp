#include "trans/t_diffs.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

TDiffs::Stats TDiffs::analyze(const std::string &f, const Options &options)
{
    TDiffs::Stats stats;
    const auto &r = Standard::instance().r_trans;

    auto c = std::map<std::string, Counts>(); //TODO (options.level == Gene ? Analyzer::baseHist() : Analyzer::seqHist());

    options.info("Parsing input file");

    // This is needed to determine any undetected sequin
    std::set<SequinID> ids;
    
    auto g = [&](const GeneID &id, double fpkm_1, double fpkm_2)
    {
        // The known and observed fold-change
        Fold known = NAN;
        
        // It's NAN if the sequin defined in reference but not in mixture
        Fold measured = NAN;
        
        const auto *g = r.findGene(id);
        
        if (g)
        {
            // Calculate the known fold-change between B and A
            known = (g->abund(MixB) / g->abund(MixA));
        }

        if (g && !isnan(fpkm_1) && !isnan(fpkm_2) && fpkm_1 && fpkm_2)
        {
            c[id]++;

            // Measured fold-change between the two mixtures
            measured = fpkm_2 / fpkm_1;
        }

        ids.insert(id);        
        stats.add(id, !isnan(known) ? log2(known) : NAN, !isnan(measured) ? log2(measured) : NAN);
    };

    ParserCDiffs::parse(f, [&](const TrackingDiffs &t, const ParserProgress &)
    {
        static const auto &id = Standard::instance().id;

        if (t.chromID != id)
        {
            stats.n_hg38++;
            return;
        }
        
        stats.n_chrT++;
        
        // The known and observed fold-change
        Fold known = NAN;

        // It's NAN if the sequin defined in reference but not in mixture
        Fold measured = NAN;

        switch (options.level)
        {
            case Gene:
            {
                if (t.status != NoTest)
                {
                    g(t.geneID, t.fpkm_1, t.fpkm_2);
                }
                else
                {
                    g(t.geneID, NAN, NAN);
                }
                
                break;
            }

            case Isoform:
            {
                const auto *seq = r.seq(t.testID);

                if (seq)
                {
                    // Known fold-change between the two mixtures
                    known = seq->abund(MixB) / seq->abund(MixA);
                }

                if (t.status != NoTest && t.fpkm_1 && t.fpkm_2)
                {
                    c[t.testID]++;

                    // Measured fold-change between the two mixtures
                    measured = t.fpkm_2 / t.fpkm_1;
                }

                stats.add(t.testID, !isnan(known) ? log2(known) : NAN, !isnan(measured) ? log2(measured) : NAN);
                break;
            }
        }
    });

    assert(!ids.empty());
    
    /*
     * Find out the undetected sequins
     */
    
//    if (options.level == Gene)
//    {
//        for (const auto &i : seq->abund)
//        {
//            const auto &id = i.first;
//            
//            // Not found in the experiment?
//            if (!ids.count(id))
//            {
//                stats.add(i.first, s.bases_2.at(id).abund() / s.bases_1.at(id).abund(), NAN);
//            }
//        }
//    }
    
    //assert(!c.empty() && !stats.x.empty());
    //assert(!stats.x.empty() && stats.x.size() == stats.y.size());

    //stats.s = Expression::analyze(c, s.r_gene(options.rMix));

    options.info("Generating statistics");
    AnalyzeReporter::linear(stats, "TransDifferent", "FPKM", options.writer);

    return stats;
}