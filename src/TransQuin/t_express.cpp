#include <stdexcept>
#include "writers/r_writer.hpp"
#include "TransQuin/t_express.hpp"
#include "parsers/parser_kallisto.hpp"
#include "parsers/parser_cufflink.hpp"
#include "parsers/parser_stringtie.hpp"

using namespace Anaquin;

typedef TExpress::Metrics  Metrics;
typedef TExpress::Software Software;

template <typename T> void update(TExpress::Stats &stats, const T &t, const TExpress::Options &o)
{
    if (t.cID != ChrT)
    {
        stats.n_geno++;
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
                        stats.add(t.id, m->concent(Mix_1), t.abund);
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
                        stats.add(t.id, m->concent(Mix_1), t.abund);
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
            case Software::Kallisto:
            {
                struct TempData : public ParserKallisto::Data
                {
                    ChrID cID = ChrT;

                    // Dummy locus, it'll never be matched
                    Locus l;
                };
                
                ParserKallisto::parse(file, [&](const ParserKallisto::Data &data, const ParserProgress &)
                {
                    TempData tmp;
                    
                    tmp.id = data.id;
                    tmp.abund = data.abund;
                    
                    update(stats, tmp, o);
                });

                break;
            }

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
    
    TExpress::generateSummary("TransExpress_summary.stats", files, stats, o, units);
    
    /*
     * 2. Generating detailed statistics for the sequins
     */
    
    TExpress::generateCSV("TransExpress_quins.stats", stats, o);
    
    /*
     * 3. Generating abundance vs abundance (single or multiple samples)
     */
    
    TExpress::generateRAbund("TransExpress_express.R", "TransExpress_quins.stats", stats, o);
    
    /*
     * 4. Generating major plot (but only if we have isoforms...)
     */

    if (stats.size() >= 2 && o.metrs == TExpress::Metrics::Isoform)
    {
        TExpress::generateRSplice("TransExpress_splice.R", "TransExpress_quins.stats", o);
    }
}