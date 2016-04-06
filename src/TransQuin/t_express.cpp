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

static Scripts writeCSV(const std::vector<TExpress::Stats> &stats)
{
    std::set<SequinID> seqs;

    // This is the data structure that will be convenient
    std::map<unsigned, std::map<SequinID, Concent>> data;
    
    // Expected concentration
    std::map<SequinID, Concent> expect;
    
    std::stringstream ss;
    ss << "Sequin\tEAbund";

    for (auto i = 0; i < stats.size(); i++)
    {
        ss << ((boost::format("\tA%1%") % (i+1)).str());
        
        for (const auto &j : stats[i].data)
        {
            seqs.insert(j.first);
            expect[j.first]  = j.second.x;
            data[i][j.first] = j.second.y;
        }
    }
    
    ss << "\n";

    for (const auto &seq : seqs)
    {
        ss << ((boost::format("%1%\t%2%") % seq % expect.at(seq)).str());
        
        for (auto i = 0; i < stats.size(); i++)
        {
            if (data[i].count(seq))
            {
                ss << "\t" << data[i][seq];
            }
            else
            {
                ss << "\tNA";
            }
        }
        
        ss << "\n";
    }
    
    return ss.str();
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
    
    std::vector<SequinHist>   sHists;
    std::vector<LinearStats>  lStats;
    std::vector<MappingStats> mStats;
    
    for (auto i = 0; i < files.size(); i++)
    {
        mStats.push_back(stats[i]);
        sHists.push_back(stats[i].hist);
        lStats.push_back(stats[i].data);
    }

    /*
     * Generating summary statistics (single sample or replicates)
     */
    
    o.info("Generating TransExpress_summary.stats");
    
    if (files.size() == 1)
    {
        o.writer->open("TransExpress_summary.stats");
        o.writer->write(generateSummary(stats[0], files[0], units, o));
        o.writer->close();
    }
    else
    {
        o.info("Generating TransReplicate_summary.stats");
        o.writer->open("TransReplicate_summary.stats");
        o.writer->write(StatsWriter::inflectSummary(o.rChrT, o.rEndo, files, sHists, mStats, lStats, units));
        o.writer->close();
    }
    
    /*
     * Generating detailed statistics for the sequins
     */
    
    o.info("Generating TransReplicate_quin.csv");
    o.writer->open("TransReplicate_quin.csv");
    
    if (files.size() == 1)
    {
        o.writer->write(StatsWriter::writeCSV(stats[0].data, "EAbund", "MAbund"));
    }
    else
    {
        o.writer->write(writeCSV(stats));
        o.writer->close();
    }
    
    o.writer->close();
    
    /*
     * Generating abundance vs abundance (single sample or replicates)
     */
    
    o.info("Generating TransReplicate_abundAbund.R");
    o.writer->open("TransExpress_abundAbund.R");
    
    if (files.size() == 1)
    {
        o.writer->write(RWriter::createScript("TransExpress_quins.csv", PlotTAbundAbund()));
    }
    else
    {
        o.writer->write(RWriter::createScript("TransReplicate_quin.csv", PlotRAbundAbund()));
    }

    o.writer->close();
    
    /*
     * Generating major plot for all samples (but only if we have the isoforms...)
     */
    
    if (files.size() >= 2 && o.metrs == TExpress::Metrics::Isoform)
    {
        o.writer->open("TransReplicate_major.R");
        o.writer->write(RWriter::createScript("TransReplicate_quin.csv", PlotMajor()));
        o.writer->close();
    }
}