#include "trans/t_diffs.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

TDiffs::Stats TDiffs::analyze(const std::string &f, const Options &options)
{
    TDiffs::Stats stats;
    const auto &s = Standard::instance();

    auto c = (options.level == Gene ? Analyzer::baseHist() : Analyzer::seqHist());

    options.info("Parsing input file");

    // This is needed to determine any undetected sequin
    std::set<SequinID> ids;
    
    auto g = [&](const GeneID &id, double fpkm_1, double fpkm_2)
    {
        // The known and observed fold-change
        Fold known = NAN;
        
        // It's NAN if the sequin defined in reference but not in mixture
        Fold measured = NAN;
        
        if (s.bases_1.count(id))
        {
            // Calculate the known fold-change between B and A
            known = (s.bases_2.at(id).abund() / s.bases_1.at(id).abund());
        }

        if (s.bases_1.count(id) && !isnan(fpkm_1) && !isnan(fpkm_2) && fpkm_1 && fpkm_2)
        {
            c[id]++;

            // Measured fold-change between the two mixtures
            measured = fpkm_2 / fpkm_1;
        }

        ids.insert(id);
        stats.z.push_back(id);
        stats.x.push_back(!isnan(known)    ? log2(known)    : NAN);
        stats.y.push_back(!isnan(measured) ? log2(measured) : NAN);
    };
    
    ParserCDiffs::parse(f, [&](const TrackingDiffs &t, const ParserProgress &)
    {
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
                if (s.seqs_1.count(t.testID))
                {
                    // Known fold-change between the two mixtures
                    known = s.seqs_2.at(t.testID).abund() / s.seqs_1.at(t.testID).abund();
                }
                
                if (t.status != NoTest && t.fpkm_1 && t.fpkm_2)
                {
                    c[t.testID]++;

                    // Measured fold-change between the two mixtures
                    measured = t.fpkm_2 / t.fpkm_1;
                }

                stats.z.push_back(t.testID);
                stats.x.push_back(!isnan(known)    ? log2(known)    : NAN);
                stats.y.push_back(!isnan(measured) ? log2(measured) : NAN);
                
                break;
            }
        }
    });

    assert(!ids.empty());
    
    /*
     * Find out the undetected sequins
     */
    
    if (options.level == Gene)
    {
        for (const auto &i : s.bases_1)
        {
            const auto &id = i.first;
            
            // Not found in the experiment?
            if (!ids.count(id))
            {
                stats.z.push_back(i.first);
                stats.x.push_back(s.bases_2.at(id).abund() / s.bases_1.at(id).abund());
                stats.y.push_back(NAN);
            }
        }
    }
    
    assert(!c.empty() && !stats.x.empty());
    assert(!stats.x.empty() && stats.x.size() == stats.y.size());

    //stats.s = Expression::analyze(c, s.r_gene(options.rMix));

    options.info("Generating statistics");
    AnalyzeReporter::linear(stats, "TransDifferent", "FPKM", options.writer);

    return stats;
}