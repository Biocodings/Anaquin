#include "trans/t_diffs.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

TDiffs::Stats TDiffs::analyze(const std::string &file, const Options &o)
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

        stats.add(id, !isnan(known) ? log2(known) : NAN, !isnan(measured) ? log2(measured) : NAN);
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
                    known = seq->abund(Mix_2) / seq->abund(Mix_1);
                }

                if (t.status != NoTest && t.fpkm_1 && t.fpkm_2)
                {
                    stats.h.at(t.testID)++;

                    // Measured fold-change between the two mixtures
                    measured = t.fpkm_2 / t.fpkm_1;
                }

                stats.add(t.testID, !isnan(known) ? log2(known) : NAN, !isnan(measured) ? log2(measured) : NAN);
                break;
            }
        }
    });

    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Genome: %2% reads\n"
                         "   Query: %3% reads\n"
                         "   Reference: %4% sequins\n\n"
                         "   Detected: %5%\n\n"
                         "   Correlation:\t%6%\n"
                         "   Slope:\t%7%\n"
                         "   R2:\t%8%\n"
                         "   Adjusted R2:\t%9%\n"
                         "   F-statistic:\t%10%\n"
                         "   P-value:\t%11%\n"
                         "   SSM: %12%, DF: %13%\n"
                         "   SSE: %14%, DF: %15%\n"
                         "   SST: %16%, DF: %17%\n";
    
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
                                            % lm.ar2
                                            % lm.f
                                            % lm.p
                                            % lm.ssm
                                            % lm.ssm_df
                                            % lm.sse
                                            % lm.sse_df
                                            % lm.sst
                                            % lm.sst_df).str());
    o.writer->close();
    
    /*
     * Generate R scripts
     */
    
    AnalyzeReporter::scatter(stats, "TransDiff", "", o.writer);
    
    return stats;
}