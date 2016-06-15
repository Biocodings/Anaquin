#include <stdexcept>
#include "writers/r_writer.hpp"
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"

using namespace Anaquin;

typedef TExpress::Metrics  Metrics;

template <typename T> void update(TExpress::Stats &stats, const T &x, const TExpress::Options &o)
{
    if (Standard::isSynthetic(x.cID))
    {
        stats.n_syn++;
        const auto &r = Standard::instance().r_trans;
        
        switch (o.metrs)
        {
            case Metrics::Isoform:
            {
                const TransData *m = nullptr;
                
                // Try to match by name if possible
                m = r.match(x.id);
                
                if (!m)
                {
                    // Try to match by locus (de-novo assembly)
                    m = r.match(x.l, Overlap);
                }
                
                if (m)
                {
                    stats.hist.at(m->id)++;
                    
                    if (x.abund)
                    {
                        stats.add(m->id, m->concent(Mix_1), x.abund);
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown isoform.") % x.id).str());
                }
                
                break;
            }
                
            case Metrics::Gene:
            {
                const TransRef::GeneData *m = nullptr;
                
                // Try to match by name if possible
                m = r.findGene(x.cID, x.id);
                
                if (!m)
                {
                    // Try to match by locus (de-novo assembly)
                    m = r.findGene(x.cID, x.l, Contains);
                }
                
                if (m)
                {
                    stats.hist.at(m->id)++;
                    
                    if (x.abund)
                    {
                        stats.add(m->id, m->concent(Mix_1), x.abund);
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found.") % x.id).str());
                }
                
                break;
            }
        }
    }
    else
    {
        stats.n_gen++;
    }
}

template <typename Functor> TExpress::Stats calculate(const TExpress::Options &o, Functor f)
{
    TExpress::Stats stats;
    
    const auto &r = Standard::instance().r_trans;
    
    switch (o.metrs)
    {
        case Metrics::Isoform: { stats.hist = r.hist();         break; }
        case Metrics::Gene:    { stats.hist = r.geneHist(ChrT); break; }
    }
    
    assert(!stats.hist.empty());
    
    f(stats);
    
    if (stats.empty())
    {
        throw std::runtime_error("Failed to find anything for the synthetic chromosome");
    }
    
    switch (o.metrs)
    {
        case Metrics::Isoform: { stats.limit = r.absolute(stats.hist);     break; }
        case Metrics::Gene:    { stats.limit = r.absoluteGene(stats.hist); break; }
    }
    
    return stats;
}

TExpress::Stats TExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);
    
    return calculate(o, [&](TExpress::Stats &stats)
    {
        switch (o.inputs)
        {
            case Inputs::Text:
            {
                ParserExpress::parse(Reader(file), [&](const ParserExpress::Data &x, const ParserProgress &)
                {
                    update(stats, x, o);
                });
                
                break;
            }
                
            case Inputs::GTF:
            {
                ParserExpress::Data t;
                
                ParserGTF::parse(file, [&](const ParserGTF::Data &x, const std::string &, const ParserProgress &)
                {
                    t.l     = x.l;
                    t.cID   = x.cID;
                    t.id    = o.metrs == Metrics::Gene ? x.gID : x.tID;
                    t.abund = x.fpkm;

                    update(stats, t, o);
                });

                break;
            }
        }
    });
}

void TExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto m = std::map<TExpress::Metrics, std::string>
    {
        { TExpress::Metrics::Gene,    "genes"    },
        { TExpress::Metrics::Isoform, "isoforms" },
    };
    
    const auto units = m.at(o.metrs);
    const auto stats = analyze(files, o);

    /*
     * Generating RnaExpress_summary.stats
     */
    
    TExpress::generateSummary("RnaExpress_summary.stats", files, stats, o, units);
    
    /*
     * Generating RnaExpress_quins.csv
     */
    
    TExpress::generateCSV("RnaExpress_quins.csv", stats, o);
    
    /*
     * Generating RnaExpress_express.R
     */
    
    TExpress::generateR("RnaExpress_express.R", "RnaExpress_quins.csv", stats, o);
}