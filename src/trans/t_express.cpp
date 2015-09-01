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
    const auto &s = Standard::instance();

    // Detect whether it's a file of isoform by the name of file
    const bool isoform = options.level == Isoform;

    options.info("Parsing input file");
    options.logInfo(isoform ? "Isoform tracking" : "Gene tracking");

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);

        options.logInfo((boost::format("%1%: %2%") % p.i % t.trackID).str());

        if (isoform)
        {
            // Try to match by names if possible
            const auto *r = s.seqs_1.count(t.trackID) ? &(s.seqs_1.at(t.trackID)) : nullptr;

            if (!r)
            {
                const auto found = std::find_if(s.seqs_1.begin(), s.seqs_1.end(),
                                                [&](const std::pair<SequinID, Sequin> &p)
                                                {
                                                    return p.second.l.contains(t.l);
                                                });

                if (found != s.seqs_1.end())
                {
                    r = &(found->second);
                }
            }

            if (!r)
            {
                options.logWarn((boost::format("%1% not found. Unknown isoform.") % t.trackID).str());
            }
            else
            {
                stats.c[t.trackID]++;
                
                if (t.fpkm)
                {
                    stats.add(t.trackID, log2(r->abund() / r->length), log2(fpkm));
                }
            }
        }
        else
        {
            // Try to match by names if possible
            const auto *r = s.bases_1.count(t.geneID) ? &(s.bases_1.at(t.geneID)) : nullptr;

            if (!r)
            {
                const auto found = std::find_if(s.seqs_1.begin(), s.seqs_1.end(),
                                        [&](const std::pair<SequinID, Sequin> &p)
                                        {
                                            return p.second.l.contains(t.l);
                                        });

                if (found != s.seqs_1.end())
                {
                    r = &(s.bases_1.at(found->second.baseID));
                }
            }

            if (!r)
            {
                options.logWarn((boost::format("%1% not found. Unknown gene.") % t.trackID).str());
            }
            else
            {
                stats.c[t.geneID]++;

                if (t.fpkm)
                {
                    stats.add(t.trackID, log2(r->abund() / r->length()), log2(fpkm));
                }
            }
        }
    });
    
    if (isoform)
    {
        stats.s = Expression::analyze(stats.c, s.seqs_1);
    }
    else
    {
        stats.s = Expression::analyze(stats.c, s.bases_1);
    }

    {
        options.info("Generating statistics");
        AnalyzeReporter::linear(stats, "TransExpression", "FPKM", options.writer);
    }
        
    return stats;
}