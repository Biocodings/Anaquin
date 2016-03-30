#include <stdexcept>
#include "data/experiment.hpp"
#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_kallisto.hpp"
#include "parsers/parser_tracking.hpp"
#include "parsers/parser_stringtie.hpp"
#include <ss/regression/segmented.hpp>

using namespace Anaquin;

typedef TExpress::Metrics  Metrics;
typedef TExpress::Software Software;

struct InternalKallistoData : public ParserKallisto::Data
{
    ChrID cID = ChrT;
    
    // Dummy value...
    Locus l;
};

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
                        stats.data[t.cID].add(t.id, m->abund(Mix_1), t.abund);
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
                        stats.data[t.cID].add(t.id, m->abund(Mix_1), t.abund);
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

template <typename Functor> TExpress::Stats calculate(const SampleName &name, const TExpress::Options &o, Functor f)
{
    TExpress::Stats stats;
    
    // Eg: A1
    stats.name = name;
    
    const auto &r   = Standard::instance().r_trans;
    const auto cIDs = r.ChrIDs();
    
    stats.data[ChrT];
    stats.data[Endo];
    
    switch (o.metrs)
    {
        case Metrics::Isoform: { stats.hist = r.hist();         break; }
        case Metrics::Gene:    { stats.hist = r.geneHist(ChrT); break; }
    }
    
    f(stats);
    
    if (stats.data.at(ChrT).empty())
    {
        throw std::runtime_error("Failed to find anything from the synthetic chromosome");
    }
    
    switch (o.metrs)
    {
        case Metrics::Isoform: { stats.limit = r.absolute(stats.hist);  break; }
        case Metrics::Gene:    { stats.limit = r.limitGene(stats.hist); break; }
    }
    
    return stats;
}

TExpress::Stats TExpress::analyze(const std::vector<Expression> &exps, const Options &o)
{
    return calculate("...", o, [&](TExpress::Stats &stats)
    {
        for (const auto &i : exps)
        {
            update(stats, i, o);
        }
    });
}

TExpress::Stats TExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);

    return calculate(o.exp->fileToSample(file), o, [&](TExpress::Stats &stats)
    {
        switch (o.soft)
        {
            case Software::Kallisto:
            {
                ParserKallisto::parse(file, [&](const ParserKallisto::Data &data, const ParserProgress &p)
                {
                    InternalKallistoData x;
                    
                    x.id    = data.id;
                    x.abund = data.abund;

                    update(stats, x, o);
                });

                break;
            }

            case Software::Cufflinks:
            {
                ParserTracking::parse(file, [&](const ParserTracking::Data &data, const ParserProgress &p)
                {
                    update(stats, data, o);
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

static void writeSummary(const TExpress::Stats &stats,
                         const FileName        &file,
                         const std::string     &name,
                         const Units           &units,
                         const TExpress::Options &o)
{
    o.writer->create(name);
    o.writer->open(name + "/TransExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                std::vector<FileName>     { file  },
                                                std::vector<SequinHist>   { stats.hist },
                                                std::vector<MappingStats> { stats },
                                                std::vector<LinearStats>  { stats.data.at(ChrT) },
                                                units));
    o.writer->close();
}

template <typename Stats, typename Options> Scripts writeSampleCSV(const std::vector<SequinID> &ids, const Stats &stats, const Options &o)
{
    assert(!ids.empty());
    
    std::stringstream ss;
    
    /*
     * 1: Generating the sample names
     */
    
    const auto &names = o.exp->names();
    
    for (const auto &name : names)
    {
        ss << ("," + name);
    }
    
    ss << "\n";
    
    /*
     * 2: Generating for the features
     */

    // Expected for each sequin across samples (should be identical)
    std::map<SequinID, std::map<SampleName, double>> expected;

    // Measured for each sequin across samples
    std::map<SequinID, std::map<SampleName, double>> measured;

    for (const auto &i : stats)
    {
        for (const auto &j : i.data.at(ChrT))
        {
            // Eg: R1_1_1
            const auto &id = j.first;
            
            // Eg: A1
            const auto name = i.name;

            expected[id][i.name] = j.second.x;
            measured[id][i.name] = j.second.y;
        }
    }
    
    for (const auto &id : ids)
    {
        ss << id;
        
        if (!expected.count(id))
        {
            ss << ",NA";
            
            for (auto i = 0; i < names.size(); i++)
            {
                ss << ",NA";
            }
        }
        else
        {
            /*
             * Generating expected concentration (any replicate will give the identical concentration)
             */
            
            ss << "," << expected[id][names.front()];
            
            /*
             * Generating measured concentration for all samples
             */
            
            for (const auto &name : names)
            {
                if (!measured[id].count(name))
                {
                    ss << ",NA";
                }
                else
                {
                    ss << "," << measured[id].at(name);
                }
            }
        }
        
        ss << "\n";
    }
    
    return ss.str();
}

static void writeFPKM(const FileName &file, const std::vector<TExpress::Stats> &stats, const TExpress::Options &o)
{
    o.writer->open(file);
    o.writer->write("expected", false);

    /*
     * Generating sample headers
     */
    
    const auto &names = o.exp->names();
    
    assert(names.size() == stats.size());
    
    for (const auto &name : names)
    {
        o.writer->write("," + name, false);
    }
    
    o.writer->write("\n", false);

    // Number of samples
    const auto n = names.size();
    
    auto f = [&](const ChrID &id)
    {
        for (const auto &i : stats[0].data.at(id))
        {
            o.writer->write(i.first, false);
            
            /*
             * Generating expected expression
             */

            o.writer->write("," + std::to_string(stats[0].data.at(id).at(i.first).x), false);

            /*
             * Generating measured expression
             */

            // For all the samples...
            for (auto j = 0; j < n; j++)
            {
                if (stats[j].data.at(id).count(i.first))
                {
                    o.writer->write("," + std::to_string(stats[j].data.at(id).at(i.first).y), false);
                }
                else
                {
                    o.writer->write(",0", false);
                }
            }
            
            o.writer->write("\n", false);
        }
    };
    
    // Writing features for synthetic chromosome
    f(ChrT);

    // Writing features for endogenous
    f(Endo);
    
    o.writer->close();
}

static void writeScatter(const TExpress::Stats   &stats,
                         const FileName          &file,
                         const std::string       &name,
                         const std::string       &units,
                         const TExpress::Options &o)
{
    o.writer->create(name);
    o.writer->open(name + "/TransExpress_scatter.R");
    o.writer->write(RWriter::scatter(stats, ChrT, "????",
                                     "TransExpress",
                                     "Expected concentration (attomol/ul)",
                                     "Measured coverage (FPKM)",
                                     "Expected concentration (log2 attomol/ul)",
                                     "Measured coverage (log2 FPKM)"));
    o.writer->close();
}

void TExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TExpress::analyze(files, o);
    
    const auto m = std::map<TExpress::Metrics, std::string>
    {
        { TExpress::Metrics::Gene,    "gene"    },
        { TExpress::Metrics::Isoform, "isoform" },
    };
    
    const auto units = m.at(o.metrs);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics for each sample
     */
    
    std::vector<TExpress::Stats::Data> data;
    std::vector<MappingStats> data_;
    std::vector<SequinHist> hists;
    
    for (auto i = 0; i < files.size(); i++)
    {
        // Generating summary statistics for the sample
        writeSummary(stats[i], files[i], o.exp->names().at(i), units, o);
        
        o.writer->open(o.exp->names().at(i) + "/TransExpress_quins.csv");
        o.writer->write(StatsWriter::writeCSV(stats[i].data.at(ChrT), "Expected concentration (attomol/ul)", "Measured abundance (attomol/ul)"));
        o.writer->close();
        
        // Generating scatter plot for the sample
        writeScatter(stats[i], files[i], o.exp->names().at(i), units, o);

        data_.push_back(stats[i]);
        hists.push_back(stats[i].hist);
        data.push_back(stats[i].data.at(ChrT));
    }
    
    /*
     * Generating CSV of expression for all samples
     */
  
//    const auto &r = Standard::instance().r_trans;
//
//    o.writer->open("ABCD");
//    
//    switch (o.lvl)
//    {
//        case Level::Gene:
//        {
//            o.writer->write(writeSampleCSV(r.geneIDs(ChrT), stats, o));
//            break;
//        }
//
//        case Level::Isoform:
//        {
//            o.writer->write(writeSampleCSV(r.seqIDs(), stats, o));
//            break;
//        }
//
//        case Level::Exon:
//        {
//            break;
//        }
//    }
//
//    o.writer->close();

    /*
     * Generating a CSV of expression for all samples
     */
    
    writeFPKM("TExpress_FPKM.csv", stats, o);

    /*
     * Generating summary statistics for all samples
     */
    
    o.writer->open("TransExpress_pooled.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT, o.rEndo, files, hists, data_, data, units));
    o.writer->close();

    /*
     * Generating scatter plot for all samples
     */
    
    o.writer->open("TransExpress_pooled.R");
    o.writer->write(RWriter::scatterPool("TExpress_FPKM.csv"));
    o.writer->close();

    /*
     * Generating spliced plot for all samples (but only if we have the isoforms...)
     */

    if (o.metrs == TExpress::Metrics::Isoform)
    {
        o.writer->open("TransExpress_Splice.R");
        o.writer->write(RWriter::createSplice("TExpress_FPKM.csv"));
        o.writer->close();
    }
}