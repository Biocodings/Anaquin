#include <stdexcept>
#include "writers/r_writer.hpp"
#include "TransQuin/t_express.hpp"
#include "parsers/parser_cufflink.hpp"
#include "parsers/parser_stringtie.hpp"

using namespace Anaquin;

extern Scripts PlotTAbundAbund();

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
                        stats.data.add(t.id, m->abund(Mix_1), t.abund);
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
                        stats.data.add(t.id, m->abund(Mix_1), t.abund);
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
    
    if (stats.data.empty())
    {
        throw std::runtime_error("Failed to find anything for the synthetic chromosome");
    }
    
    switch (o.metrs)
    {
        case Metrics::Isoform: { stats.data.limit = r.absolute(stats.hist);     break; }
        case Metrics::Gene:    { stats.data.limit = r.absoluteGene(stats.hist); break; }
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
                ParserCufflink::parse(file, [&](const ParserCufflink::Data &data, const ParserProgress &p)
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

static Scripts generateSummary(const TExpress::Stats &stats, const FileName &file, const Units &units, const TExpress::Options &o)
{
    return StatsWriter::inflectSummary(o.rChrT,
                                       o.rEndo,
                                       std::vector<FileName>     { file  },
                                       std::vector<SequinHist>   { stats.hist },
                                       std::vector<MappingStats> { stats },
                                       std::vector<LinearStats>  { stats.data },
                                       units);
}

void TExpress::report(const FileName &file, const Options &o)
{
    const auto m = std::map<TExpress::Metrics, std::string>
    {
        { TExpress::Metrics::Gene,    "genes"    },
        { TExpress::Metrics::Isoform, "isoforms" },
    };
    
    const auto stats = TExpress::analyze(file, o);
    const auto units = m.at(o.metrs);
    
    /*
     * Generating summary statistics (all together if multiple samples)
     */

    o.info("Generating TransExpress_summary.stats");
    o.writer->open("TransExpress_summary.stats");
    o.writer->write(generateSummary(stats, file, units, o));
    o.writer->close();
    
    /*
     * Generating detailed statistics for the sequins
     */
    
    o.info("Generated TransExpress_quins.csv");
    o.writer->open("TransExpress_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats.data, "EAbund", "MAbund"));
    o.writer->close();
    
    /*
     * Generating for AbundAbund
     */
    
    o.info("Generated TranExpress_abundAbund.R");
    o.writer->open("TransExpress_abundAbund.R");
    o.writer->write(RWriter::createScript("TransExpress_quins.csv", PlotTAbundAbund()));
    o.writer->close();
}