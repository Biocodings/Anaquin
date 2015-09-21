#include "trans/t_diffs.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

TDiffs::Stats TDiffs::report(const std::string &file, const Options &o)
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

    o.info("Generating statistics");

    const auto units = isoform ? "isoforms" : "genes";

    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Genome: %2% %17%\n"
                         "   Query: %3% %17%\n"
                         "   Reference: %4% %17%\n\n"
                         "   Detected: %5% %17%\n\n"
                         "   Correlation:\t%6%\n"
                         "   Slope:\t%7%\n"
                         "   R2:\t%8%\n"
                         "   F-statistic:\t%9%\n"
                         "   P-value:\t%10%\n"
                         "   SSM: %11%, DF: %12%\n"
                         "   SSE: %13%, DF: %14%\n"
                         "   SST: %15%, DF: %16%\n";
    
    const auto lm = stats.linear();

    o.writer->open("TransDiff_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.n_hg38
                                            % stats.n_chrT
                                            % stats.h.size()
                                            % countHist(stats.h)
                                            % lm.r
                                            % lm.m
                                            % lm.r2
                                            % lm.f
                                            % lm.p
                                            % lm.ssm
                                            % lm.ssm_df
                                            % lm.sse
                                            % lm.sse_df
                                            % lm.sst
                                            % lm.sst_df
                                            % units).str());
    o.writer->close();
    
    /*
     * Generate an R script
     */
    
    AnalyzeReporter::scatter(stats, "TransDiff", "Expected log-fold of mixture A and B", "Measured log-fold of mixture A and B", o.writer);
    
    return stats;
}