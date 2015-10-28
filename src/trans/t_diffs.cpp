#include "trans/t_diffs.hpp"
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

TDiffs::Stats TDiffs::report(const FileName &file, const Options &o)
{
    TDiffs::Stats stats;
    const auto &r = Standard::instance().r_trans;

    const bool isoform = o.level == Isoform;
    o.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    // Construct for a histogram at the appropriate level
    stats.h = isoform ? r.hist() : r.histGene();
    
    o.info("Parsing tracking file");

    /*
     * Differential expression at the gene level
     */
    
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
            known = (g->abund(Mix_2) / g->abund(Mix_1));
        }

        if (g && !isnan(fpkm_1) && !isnan(fpkm_2) && fpkm_1 && fpkm_2)
        {
            stats.h.at(id)++;

            // Measured fold-change between the two mixtures
            measured = fpkm_2 / fpkm_1;
        }

        stats.add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
    };

    ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
    {
        static const auto &id = Standard::chrT;

        if (t.chromID != id)
        {
            stats.n_expT++;
            return;
        }
        
        stats.n_chrT++;
        
        // The known and observed fold-change
        Fold known = NAN;

        // It's NAN if the sequin defined in reference but not in mixture
        Fold measured = NAN;

        switch (o.level)
        {
            case Gene:
            {
                if (t.status != NoTest && stats.h.count(t.geneID))
                {
                    g(t.geneID, t.fpkm_1, t.fpkm_2);
                }
                
                break;
            }

            case Isoform:
            {
                if (t.status == NoTest || !stats.h.count(t.testID))
                {
                    return;
                }
                
                const auto *seq = r.match(t.testID);

                if (seq)
                {
                    // Known fold-change between the two mixtures
                    known = seq->abund(Mix_2) / seq->abund(Mix_1);
                }

                if (t.status != NoTest && t.fpkm_1 && t.fpkm_2)
                {
                    stats.h.at(t.testID)++;

                    // Measured fold-change between the two mixtures
                    measured = t.fpkm_2 / t.fpkm_1;
                }

                stats.add(t.testID, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);

                break;
            }
        }
    });

    o.info("Calculating limit of sensitivity");

    stats.ss = isoform ? r.limit(stats.h) : r.limitGene(stats.h);

    const auto units = isoform ? "isoforms" : "genes";

    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    AnalyzeReporter::linear("TransDifferent_summary.stats", file, stats, units, o.writer);

    /*
     * Generating Bioconductor
     */
    
    o.info("Generating Bioconductor");
    AnalyzeReporter::scatter(stats, "", "TransDifferent", "Expected fold change of mixture A and B", "Measured fold change of mixture A and B", "Expected log2 fold change of mixture A and B", "Expected log2 fold change of mixture A and B", o.writer);
    
    return stats;
}