#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_tracking.hpp"

using namespace SS;
using namespace Anaquin;

TExpress::Stats TExpress::analyze(const std::string &file, const Options &o)
{
    TExpress::Stats stats;
    const auto &r = Standard::instance().r_trans;

    const bool isoform = o.level == Isoform;
    o.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    // Construct for a histogram at the appropriate level
    stats.h = isoform ? r.hist() : r.histGene();

    o.info("Parsing input file");

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
            m = r.seq(t.geneID);

            if (!m)
            {
                // Try to match by locus (de-novo assembly)
                m = r.seq(t.l);
            }

            if (!m)
            {
                o.logWarn((boost::format("%1% not found. Unknown isoform.") % t.trackID).str());
            }
            else
            {
                stats.h.at(t.geneID)++;

                if (t.fpkm)
                {
                    stats.add(t.trackID, log2(m->abund(Mix_1)), log2(fpkm));
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
                m = r.findGene(t.l, TransRef::Exact);
            }

            if (!m)
            {
                o.logWarn((boost::format("%1% not found. Unknown gene.") % t.trackID).str());
            }
            else
            {
                stats.h.at(t.geneID)++;

                if (t.fpkm)
                {
                    stats.add(t.trackID, log2(m->abund(o.mix)), log2(fpkm));
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
                         "   SST: %16%, DF: %17%\n"
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
    
    AnalyzeReporter::scatter(stats, "TransExpress", "", o.writer);

    return stats;
}