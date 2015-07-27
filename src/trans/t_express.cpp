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

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);

        options.logInfo((boost::format("%1%: %2%") % p.i % t.trackID).str());

        if (isoform)
        {
            // Try to match by names if possible
            const auto *r = s.r_seqs_A.count(t.trackID) ? &(s.r_seqs_A.at(t.trackID)) : nullptr;

            if (!r)
            {
                const auto found = std::find_if(s.r_seqs_A.begin(), s.r_seqs_A.end(),
                                                [&](const std::pair<SequinID, Sequin> &p)
                                                {
                                                    return p.second.l.contains(t.l);
                                                });

                if (found != s.r_seqs_A.end())
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
                    stats.x.push_back(log2(r->abund() / r->length));
                    stats.y.push_back(log2(fpkm));
                    stats.z.push_back(t.trackID);
                }
            }
        }
        else
        {
            // Try to match by names if possible
            const auto *r = s.r_seqs_gA.count(t.geneID) ? &(s.r_seqs_gA.at(t.geneID)) : nullptr;

            if (!r)
            {
                const auto found = std::find_if(s.r_seqs_A.begin(), s.r_seqs_A.end(),
                                                [&](const std::pair<SequinID, Sequin> &p)
                                                {
                                                    return p.second.l.contains(t.l);
                                                });
                
                if (found != s.r_seqs_A.end())
                {
                    r = &(s.r_seqs_gA.at(found->second.baseID));
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
                    stats.x.push_back(log2(r->abund() / r->length()));
                    stats.y.push_back(log2(fpkm));
                    stats.z.push_back(t.trackID);
                }
            }
        }
    });
    
    assert(!stats.x.empty() && stats.x.size() == stats.y.size() && stats.y.size() == stats.z.size());

    if (isoform)
    {
        stats.s = Expression::analyze(stats.c, s.r_sequin(options.mix));
    }
    else
    {
        stats.s = Expression::analyze(stats.c, s.r_gene(options.mix));
    }
    
    options.info("Generating statistics");

    AnalyzeReporter::linear(stats, "TransExpression", "FPKM", options.writer);

    return stats;
}