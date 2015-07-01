#include <iostream>
#include <ss/c.hpp>
#include "rna/r_abundance.hpp"
#include "writers/r_writer.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_tmap.hpp"
#include "parsers/parser_tracking.hpp"

using namespace SS;
using namespace Spike;

static const FileName GTracking = "genes.fpkm_tracking";
static const FileName ITracking = "isoforms.fpkm_tracking";

RAbundanceStats RAbundance::analyze(const std::string &file, const Options &options)
{
    RAbundanceStats stats;
    const auto &s = Standard::instance();

    auto c = RAnalyzer::sequinCounter();
    unsigned i = 0;
    
    if (file.find(GTracking) == std::string::npos && file.find(ITracking) == std::string::npos)
    {
        throw std::invalid_argument((boost::format("Unknown file. It must be %1% or %2%")
                                     % GTracking
                                     % ITracking).str());
    }
    
    // Detect whether it's a file of isoform by the name of file
    const bool isoform = file.find(ITracking) != std::string::npos;
    
    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &)
    {
        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);

        if (isoform)
        {
            if (!s.r_seqs_A.count(t.trackID))
            {
                return;
            }
            
            i++;
            c[t.trackID]++;
            
            if (t.fpkm)
            {
                const auto &i = s.r_seqs_A.at(t.trackID);

                stats.x.push_back(log(i.abund() / i.length));
                stats.y.push_back(log(fpkm));
                stats.z.push_back(t.trackID);
            }
            
            stats.s = Expression::analyze(c, s.r_sequin(options.mix));
        }
        else
        {
            if (!s.r_seqs_gA.count(t.trackID))
            {
                return;
            }
            
            i++;
            c[t.geneID]++;
            
            const auto &m = s.r_seqs_gA.at(t.geneID);
            
            if (t.fpkm)
            {
                stats.x.push_back(log(m.abund()));
                stats.y.push_back(log(fpkm));
                stats.z.push_back(t.geneID);
            }
            
            stats.s = Expression::analyze(c, s.r_gene(options.mix));
        }
    });

    assert(!stats.x.empty() && stats.x.size() == stats.y.size() && stats.y.size() == stats.z.size());

    // Perform a linear regreession
    stats.linear();

    /*
     * Write out the statistics
     */
    
    AnalyzeReporter::report("rna_abundance.stats", "abundance.R", stats, "FPKM", c, options.writer);

    /*
     * Write out results for RNA sequins
     */

    const std::string format = "%1%\t%2%\t%3%";
    
    options.writer->open("rna_sequins.stats");
    options.writer->write((boost::format(format) % "id" % "spiked" % "measured").str());
    
    for (std::size_t i = 0; i < stats.z.size(); i++)
    {
        options.writer->write((boost::format(format) % stats.z[i] % stats.x[i] % stats.y[i]).str());
    }

    options.writer->close();
    
    return stats;
}