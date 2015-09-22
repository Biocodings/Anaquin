#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_tracking.hpp"

using namespace SS;
using namespace Anaquin;

TExpress::Stats TExpress::report(const std::string &file, const Options &o)
{
    TExpress::Stats stats;
    const auto &r = Standard::instance().r_trans;

    const bool isoform = o.level == Isoform;
    o.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    // Construct for a histogram at the appropriate level
    stats.h = isoform ? r.hist() : r.histGene();
    
    o.info("Parsing: " + file);

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
        static const auto &id = Standard::instance().id;
        
        if (t.chromID != id)
        {
            stats.n_hg38++;
            return;
        }

        stats.n_chrT++;

        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);

        if (isoform)
        {
            const TransData *m = nullptr;

            // Try to match by name if possible
            m = r.match(t.trackID);

            if (!m)
            {
                // Try to match by locus (de-novo assembly)
                m = r.match(t.l, Overlap);
            }

            if (!m)
            {
                o.logWarn((boost::format("%1% not found. Unknown isoform.") % t.trackID).str());
            }
            else
            {
                stats.h.at(m->id)++;

                if (t.fpkm)
                {
                    stats.add(t.trackID, m->abund(Mix_1), fpkm);
                }
            }
        }
        else
        {
            const TransRef::GeneData *m = nullptr;
            
            // Try to match by name if possible
            m = r.findGene(t.geneID);

            if (!m)
            {
                // Try to match by locus (de-novo assembly)
                m = r.findGene(t.l, Contains);
            }

            if (!m)
            {
                o.logWarn((boost::format("%1% not found. Unknown gene.") % t.trackID).str());
            }
            else
            {
                stats.h.at(m->id)++;

                if (t.fpkm)
                {
                    stats.add(t.trackID, m->abund(o.mix), fpkm);
                }
            }
        }
    });
    
    if (isoform)
    {
        stats.s = r.limit(stats.h);
    }
    else
    {
        stats.s = r.limitGene(stats.h);
    }
    
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
                         "   SST: %15%, DF: %16%\n"
    ;

    const auto lm = stats.linear();

    o.writer->open("TransExpress_summary.stats");
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
    
    AnalyzeReporter::scatter(stats, "", "TransExpress", "Expected concentration (attomol/ul)", "Measured coverage (FPKM)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 FPKM)", o.writer);

    return stats;
}