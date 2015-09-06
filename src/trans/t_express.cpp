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
    stats.h = isoform ? r.hist() : r.histForGene();

    o.info("Parsing input file");

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
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
                    stats.add(t.trackID, log2(m->abund(MixA)), log2(fpkm));
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
                m = r.findGene(t.l);
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
                    stats.add(t.trackID, log2(m->abund()), log2(fpkm));
                }
            }
        }
    });
    
    if (isoform)
    {
        stats.s = Expression_::calculate(stats.h, r);
    }
    else
    {
        stats.s = Expression_::calculate(stats.h, r);
    }
    
    o.info("Generating statistics");
    
    const auto detected = std::count_if(
                              stats.h.begin(), stats.h.end(), [&](const std::pair<SequinID, Counts> &i)
                              {
                                  return i.second;
                              });
    
    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Detected: %2%\n"
                         "   Reference: %3%\n\n"
                         "Correlation:     %4%\n"
                         "Slope:     %5%\n"
                         "R2:     %6%\n"
    ;

    const auto lm = stats.linear();

    o.writer->write((boost::format(summary) % file
                                            % detected
                                            % stats.h.size()
                                            % lm.r
                                            % lm.m
                                            % lm.r2).str());
    o.writer->close();
    
    return stats;
}