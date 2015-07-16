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

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);

        options.logger->write((boost::format("%1%: %2%") % p.i % t.trackID).str());
        
        if (isoform)
        {
            if (!s.r_seqs_A.count(t.trackID))
            {
                return;
            }
            
            stats.c[t.trackID]++;
            
            if (t.fpkm)
            {
                const auto &i = s.r_seqs_A.at(t.trackID);

                stats.x.push_back(log2(i.abund() / i.length));
                stats.y.push_back(log2(fpkm));
                stats.z.push_back(t.trackID);
            }

            stats.s = Expression::analyze(stats.c, s.r_sequin(options.mix));
        }
        else
        {
            if (!s.r_seqs_gA.count(t.trackID))
            {
                return;
            }
            
            stats.c[t.geneID]++;
            
            const auto &m = s.r_seqs_gA.at(t.geneID);
            const auto &r = std::find_if(s.r_seqs_A.begin(), s.r_seqs_A.end(),
                                [&](const std::pair<SequinID, Sequin> &p)
                                {
                                    return p.second.baseID == t.geneID;
                                });
            assert(r != s.r_seqs_A.end());
            
            if (t.fpkm)
            {
                stats.x.push_back(log2(m.abund()) / r->second.length);
                stats.y.push_back(log2(fpkm));
                stats.z.push_back(t.geneID);
            }
            
            stats.s = Expression::analyze(stats.c, s.r_gene(options.mix));
        }
    });

    assert(!stats.x.empty() && stats.x.size() == stats.y.size() && stats.y.size() == stats.z.size());

    options.info("Generating an R script");
    AnalyzeReporter::linear(stats, "rna_abund", "FPKM", options.writer);

    options.info("Generating statistics for sequin");
    const std::string format = "%1%\t%2%\t%3%";

    options.writer->open("rna_sequins.stats");
    options.writer->write((boost::format(format) % "id" % "expect" % "measured").str());

    for (std::size_t i = 0; i < stats.z.size(); i++)
    {
        options.writer->write((boost::format(format) % stats.z[i] % stats.x[i] % stats.y[i]).str());
    }

    options.writer->close();
    
    return stats;
}