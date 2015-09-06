#include <ss/c.hpp>
#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_tracking.hpp"

using namespace SS;
using namespace Anaquin;

TExpress::Stats TExpress::analyze(const std::string &file, const Options &options)
{
    TExpress::Stats stats;
    const auto &r = Standard::instance().r_trans;

    const bool isoform = options.level == Isoform;
    options.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    // Construct for a histogram at the appropriate level
    stats.h = isoform ? r.hist() : r.histForGene();

    options.info("Parsing input file");

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);

        if (isoform)
        {
//            // Try to match by names if possible
//            const auto *r = s.seqs_1.count(t.trackID) ? &(s.seqs_1.at(t.trackID)) : nullptr;
//
//            if (!r)
//            {
//                const auto found = std::find_if(s.seqs_1.begin(), s.seqs_1.end(),
//                                                [&](const std::pair<SequinID, Sequin> &p)
//                                                {
//                                                    return p.second.l.contains(t.l);
//                                                });
//
//                if (found != s.seqs_1.end())
//                {
//                    r = &(found->second);
//                }
//            }
//
//            if (!r)
//            {
//                options.logWarn((boost::format("%1% not found. Unknown isoform.") % t.trackID).str());
//            }
//            else
//            {
//                stats.h[t.trackID]++;
//                
//                if (t.fpkm)
//                {
//                    stats.add(t.trackID, log2(r->abund() / r->length), log2(fpkm));
//                }
//            }
        }
        else
        {
            const GeneData *m = nullptr;
            
            // Try to match by name if possible
            m = r.findGene(t.geneID);

            if (!m)
            {
                // Try to match by locus (de-novo assembly)
                m = r.findGene(t.l);
            }

            if (!m)
            {
                options.logWarn((boost::format("%1% not found. Unknown gene.") % t.trackID).str());
            }
            else
            {
                stats.h.at(t.geneID)++;

                if (t.fpkm)
                {
                    stats.add(t.trackID, log2(m->abund() / m->length()), log2(fpkm));
                }
            }
        }
    });
    
    if (isoform)
    {
       // stats.s = Expression::analyze(stats.h, s.seqs_1);
    }
    else
    {
       // stats.s = Expression::analyze(stats.h, s.bases_1);
    }

    options.info("Generating statistics");
    AnalyzeReporter::linear(stats, "TransExpression", "FPKM", options.writer);

    return stats;
}