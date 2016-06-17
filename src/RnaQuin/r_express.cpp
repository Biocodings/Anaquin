#include <stdexcept>
#include "writers/r_writer.hpp"
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"

using namespace Anaquin;

typedef RExpress::Metrics  Metrics;

template <typename T> void update(RExpress::Stats &stats, const T &x, const RExpress::Options &o)
{
    if (Standard::isSynthetic(x.cID))
    {
        stats.n_syn++;
        const auto &r = Standard::instance().r_trans;
        
        switch (o.metrs)
        {
            case Metrics::Isoform:
            {
                const auto m = r.match(x.id);
                
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
                const auto m = r.findGene(x.cID, x.id);

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

        // We'll need the information to estimate the numbers below and above the LOQ
        stats.gData[x.id].abund = x.abund;
    }
}

template <typename Functor> RExpress::Stats calculate(const RExpress::Options &o, Functor f)
{
    RExpress::Stats stats;
    
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

RExpress::Stats RExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);
    
    return calculate(o, [&](RExpress::Stats &stats)
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

void RExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto m = std::map<RExpress::Metrics, std::string>
    {
        { RExpress::Metrics::Gene,    "genes"    },
        { RExpress::Metrics::Isoform, "isoforms" },
    };
    
    const auto units = m.at(o.metrs);
    const auto stats = analyze(files, o);

    /*
     * Generating RnaExpress_summary.stats
     */
    
    RExpress::generateSummary("RnaExpress_summary.stats", files, stats, o, units);
    
    /*
     * Generating RnaExpress_quins.csv
     */
    
    RExpress::generateCSV("RnaExpress_quins.csv", stats, o);
    
    /*
     * Generating RnaExpress_express.R
     */
    
    RExpress::generateR("RnaExpress_express.R", "RnaExpress_quins.csv", stats, o);
}