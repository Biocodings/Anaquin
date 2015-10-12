#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include <ss/regression/linear.hpp>
#include "parsers/parser_tracking.hpp"

using namespace SS;
using namespace Anaquin;

TExpress::Stats TExpress::report(const FileName &file, const Options &o)
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
                    stats.add(t.trackID, m->abund(Mix_1), fpkm);
                }
            }
        }
    });
    
    stats.ss = isoform ? r.limit(stats.h) : r.limitGene(stats.h);
    
    const auto units = isoform ? "isoforms" : "genes";
    
    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    AnalyzeReporter::linear("TransExpress_summary.stats", file, stats, units, o.writer);

    /*
     * Generating an R script
     */
    
    o.info("Generating an R script");
    AnalyzeReporter::scatter(stats, "", "TransExpress", "Expected concentration (attomol/ul)", "Measured coverage (FPKM)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 FPKM)", o.writer);

    return stats;
}