#include <stdexcept>
#include "writers/r_writer.hpp"
#include "TransQuin/t_express.hpp"
#include "parsers/parser_cufflink.hpp"
#include "parsers/parser_stringtie.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotTAbundAbund();

// Defined in resources.cpp
extern Scripts PlotRAbundAbund();

// Defined in resources.cpp
extern Scripts PlotMajor();

typedef TExpress::Metrics  Metrics;
typedef TExpress::Software Software;

template <typename T> void update(TExpress::Stats &stats, const T &t, const TExpress::Options &o)
{
    if (t.cID != ChrT)
    {
        stats.n_endo++;
    }
    else
    {
        stats.n_chrT++;
    }
    
    if (t.cID == ChrT)
    {
        const auto &r = Standard::instance().r_trans;
        
        switch (o.metrs)
        {
            case Metrics::Isoform:
            {
                const TransData *m = nullptr;
                
                // Try to match by name if possible
                m = r.match(t.id);
                
                if (!m)
                {
                    // Try to match by locus (de-novo assembly)
                    m = r.match(t.l, Overlap);
                }
                
                if (!m)
                {
                    o.logWarn((boost::format("%1% not found. Unknown isoform.") % t.id).str());
                }
                else
                {
                    stats.hist.at(m->id)++;
                    
                    if (t.abund)
                    {
                        stats.add(t.id, m->abund(Mix_1), t.abund);
                    }
                }
                
                break;
            }
                
            case Metrics::Gene:
            {
                const TransRef::GeneData *m = nullptr;
                
                // Try to match by name if possible
                m = r.findGene(t.cID, t.id);
                
                if (!m)
                {
                    // Try to match by locus (de-novo assembly)
                    m = r.findGene(t.cID, t.l, Contains);
                }
                
                if (m)
                {
                    stats.hist.at(m->id)++;
                    
                    if (t.abund)
                    {
                        stats.add(t.id, m->abund(Mix_1), t.abund);
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown gene.") % t.id).str());
                }
                
                break;
            }
        }
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
        switch (o.soft)
        {
            case Software::Cufflinks:
            {
                ParserCufflink::parse(file, [&](const ParserCufflink::Data &data, const ParserProgress &)
                {
                    /*
                     * update() doesn't recognize the tID field. We'll need to replace it.
                     */
                    
                    auto tmp = data;
                    
                    if (o.metrs == Metrics::Isoform)
                    {
                        tmp.id = tmp.tID;
                    }
                    
                    update(stats, tmp, o);
                });
                
                break;
            }
                
            case Software::StringTie:
            {
                switch (o.metrs)
                {
                    case Metrics::Gene:
                    case Metrics::Isoform:
                    {
                        ParserStringTie::parseCTab(file, [&](const ParserStringTie::Data &data, const ParserProgress &)
                        {
                            update(stats, data, o);
                        });
                        
                        break;
                    }
                }
                
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
     * 1. Generating summary statistics (single or multiple samples)
     */
    
    o.info("Generating TransExpress_summary.stats");
    o.writer->open("TransExpress_summary.stats");
    
    if (files.size() == 1)
    {
        o.writer->write(singleSummary(stats[0], files[0], units, o));
    }
    else
    {
        o.writer->write(multipleSummary(files, stats, units, o));
    }

    o.writer->close();
    
    /*
     * 2. Generating detailed statistics for the sequins
     */
    
    o.info("Generating TransExpress_quins.csv");
    o.writer->open("TransExpress_quins.csv");
    
    if (files.size() == 1)
    {
        o.writer->write(StatsWriter::writeCSV(stats[0], "EAbund", "MAbund"));
    }
    else
    {
        o.writer->write(TExpress::multipleCSV(stats));
    }
    
    o.writer->close();
    
    /*
     * 3. Generating abundance vs abundance (single or multiple samples)
     */
    
    o.info("Generating TransExpress_abundAbund.R");
    o.writer->open("TransExpress_abundAbund.R");
    
    if (files.size() == 1)
    {
        o.writer->write(RWriter::createScript("TransExpress_quins.csv", PlotTAbundAbund()));
    }
    else
    {
        o.writer->write(RWriter::createScript("TransExpress_quins.csv", PlotRAbundAbund()));
    }

    o.writer->close();
    
    /*
     * 4. Generating major plot (but only if we have the isoforms...)
     */
    
    if (files.size() >= 2 && o.metrs == TExpress::Metrics::Isoform)
    {
        o.writer->open("TransExpress_major.R");
        o.writer->write(RWriter::createScript("TransExpress_quins.csv", PlotMajor()));
        o.writer->close();
    }
}