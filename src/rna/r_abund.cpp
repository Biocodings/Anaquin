#include <ss/c.hpp>
#include "rna/r_abund.hpp"
#include "writers/r_writer.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_tracking.hpp"

using namespace SS;
using namespace Anaquin;

RAbund::Stats RAbund::analyze(const std::string &file, const Options &options)
{
    RAbund::Stats stats;
    const auto &s = Standard::instance();

    // Detect whether it's a file of isoform by the name of file
    const bool isoform = options.level == Isoform;

    options.info("Parsing input file");

    // Whether the de-novo message has been shown
    bool shown = false;
    
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
                    if (!shown)
                    {
                        options.warn("Named sequins not found. De novo assembly assumed.");
                    }

                    shown = true;
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
                    if (!shown)
                    {
                        options.warn("Named sequins not found. De novo assembly assumed.");
                    }
                    
                    shown = true;
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
    
    options.info("Generating an R script");
    AnalyzeReporter::linear(stats, "rna_abund", "FPKM", options.writer);

    options.info("Generating statistics for sequin");
    const std::string format = "%1%\t%2%\t%3%";

    return stats;
}