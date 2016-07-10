#include <stdexcept>
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"

using namespace Anaquin;

typedef RExpress::Metrics Metrics;

template <typename T> void update(RExpress::Stats &stats, const T &x, const RExpress::Options &o)
{
    if (Standard::isSynthetic(x.cID))
    {
        stats.n_syn++;
        const auto &r = Standard::instance().r_trans;
        
        SequinID id;
        Concent  exp = NAN;
        Measured obs = NAN;

        switch (o.metrs)
        {
            case Metrics::Isoform:
            {
                const auto m = r.match(x.id);
                
                if (m)
                {
                    stats.hist.at(x.cID).at(m->id)++;

                    if (!isnan(x.abund) && x.abund)
                    {
                        id  = m->id;
                        exp = m->concent(Mix_1);
                        obs = x.abund;
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown sequin.") % x.id).str());
                }
                
                break;
            }

            case Metrics::Gene:
            {
                const auto m = r.findGene(x.cID, x.id);

                if (m)
                {
                    stats.hist.at(x.cID).at(x.id)++;
                    
                    if (!isnan(x.abund) && x.abund)
                    {
                        id  = x.id;
                        exp = r.concent(x.id);
                        obs = x.abund;
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown sequin gene.") % x.id).str());
                }
                
                break;
            }
        }
        
        if (!id.empty())
        {
            stats.add(id, exp, obs);
            
            if (isnan(stats.limit.abund) || exp < stats.limit.abund)
            {
                stats.limit.id = id;
                stats.limit.abund = exp;
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
        case Metrics::Isoform: { stats.hist = r.histIsof(); break; }
        case Metrics::Gene:    { stats.hist = r.histGene(); break; }
    }
    
    assert(!stats.hist.empty());
    
    f(stats);
    
    if (stats.empty())
    {
        throw std::runtime_error("Failed to find anything for the synthetic chromosome");
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
                    bool matched = false;
                    
                    switch (o.metrs)
                    {
                        case Metrics::Isoform:
                        {
                            matched = x.type == RNAFeature::Transcript;
                            break;
                        }
                            
                        case Metrics::Gene:
                        {
                            matched = x.type == RNAFeature::Gene;
                            break;
                        }
                    }
                    
                    if (matched)
                    {
                        t.l     = x.l;
                        t.cID   = x.cID;
                        t.id    = o.metrs == Metrics::Gene ? x.gID : x.tID;
                        t.abund = x.fpkm;
                        
                        update(stats, t, o);                        
                    }
                });

                break;
            }
        }
    });
}

static void writeQueries(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
}

static Scripts multipleCSV(const std::vector<RExpress::Stats> &stats)
{
    const auto &r = Standard::instance().r_trans;
    
    std::set<SequinID> seqs;
    
    // This is the data structure that will be convenient
    std::map<unsigned, std::map<SequinID, Concent>> data;
    
    // Expected concentration
    std::map<SequinID, Concent> expect;
    
    std::stringstream ss;
    ss << "ID\tLength\tExpected";
    
    for (auto i = 0; i < stats.size(); i++)
    {
        ss << ((boost::format("\tObserved%1%") % (i+1)).str());
        
        for (const auto &j : stats[i])
        {
            seqs.insert(j.first);
            expect[j.first]  = j.second.x;
            data[i][j.first] = j.second.y;
        }
    }
    
    ss << "\n";
    
    for (const auto &seq : seqs)
    {
        ss << ((boost::format("%1%\t%2%\t%3%") % seq
                % r.match(seq)->l.length()
                % expect.at(seq)).str());
        
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

static void generateCSV(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    o.info("Generating " + output);
    o.writer->open(output);
    
    if (stats.size() == 1)
    {
        const auto format = "%1%\t%2%\t%3%\t%4%";
        
        o.writer->write((boost::format(format) % "ID"
                                               % "Length"
                                               % "InputConcent"
                                               % "Observed").str());
        for (const auto &i : stats[0])
        {
            o.writer->write((boost::format(format) % i.first
                                                   % r.match(i.first)->l.length()
                                                   % i.second.x
                                                   % i.second.y).str());
        }
    }
    else
    {
        o.writer->write(multipleCSV(stats));
    }
    
    o.writer->close();
}

void RExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto m = std::map<RExpress::Metrics, std::string>
    {
        { RExpress::Metrics::Gene,    "genes"    },
        { RExpress::Metrics::Isoform, "isoforms" },
    };
    
    switch (o.metrs)
    {
        case Metrics::Gene:    { o.info("Gene Expresssion");   break; }
        case Metrics::Isoform: { o.info("Isoform Expression"); break; }
    }
    
    const auto units = m.at(o.metrs);
    const auto stats = analyze(files, o);

    for (const auto &i : stats)
    {
        o.info("Genome: " + toString(i.gData.size()));
    }

    /*
     * Generating RnaExpression_summary.stats
     */
    
    RExpress::generateSummary("RnaExpression_summary.stats", files, stats, o, units);
    
    /*
     * Generating RnaExpression_sequins.csv
     */
    
    generateCSV("RnaExpression_sequins.csv", stats, o);

    /*
     * Generating RnaExpression_queries.csv
     */
    
    writeQueries("RnaExpression_queries.csv", stats, o);
    
    /*
     * Generating RnaExpression_sequins.csv
     */
    
    RExpress::generateR("RnaExpression_express.R", "RnaExpression_sequins.csv", stats, o);
}