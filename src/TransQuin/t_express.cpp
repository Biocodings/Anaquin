#include <stdexcept>
#include "writers/r_writer.hpp"
#include "TransQuin/TransQuin.hpp"
#include "TransQuin/t_express.hpp"
#include "parsers/parser_stringtie.hpp"

//extern Scripts PlotTAbundAbund();

using namespace Anaquin;

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

TExpress::Stats TExpress::analyze(const std::vector<Expression> &exps, const Options &o)
{
    return calculate(o, [&](TExpress::Stats &stats)
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
    
    return calculate(o, [&](TExpress::Stats &stats)
    {
        switch (o.soft)
        {
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

template <typename Stats, typename Options> Scripts writeSampleCSV(const std::vector<SequinID> &ids, const Stats &stats, const Options &o)
{
    assert(!ids.empty());
    
    std::stringstream ss;
    
    /*
     * 1: Generating the sample names
     */
    
    //const auto &names = o.exp->names();
    
//    for (const auto &name : names)
    {
    //    ss << ("," + name);
    }
    
   // ss << "\n";
    
    /*
     * 2: Generating for the features
     */
    
    // Expected for each sequin across samples (should be identical)
    std::map<SequinID, std::map<SampleName, double>> expected;
    
    // Measured for each sequin across samples
    std::map<SequinID, std::map<SampleName, double>> measured;
    
    for (const auto &i : stats)
    {
        for (const auto &j : i.data)
        {
            // Eg: R1_1_1
            const auto &id = j.first;
            
            // Eg: A1
            //const auto name = i.name;
            
            //expected[id][i.name] = j.second.x;
            //measured[id][i.name] = j.second.y;
        }
    }
    
    for (const auto &id : ids)
    {
        ss << id;
        
        if (!expected.count(id))
        {
            ss << ",NA";
            
           // for (auto i = 0; i < names.size(); i++)
            {
                ss << ",NA";
            }
        }
        else
        {
            /*
             * Generating expected concentration (any replicate will give the identical concentration)
             */
            
           // ss << "," << expected[id][names.front()];
            
            /*
             * Generating measured concentration for all samples
             */
            
           // for (const auto &name : names)
            {
           //     if (!measured[id].count(name))
                {
                    ss << ",NA";
                }
            //    else
                {
                //    ss << "," << measured[id].at(name);
                }
            }
        }
        
        ss << "\n";
    }
    
    return ss.str();
}

static void writeScatter(const TExpress::Stats   &stats,
                         const FileName          &file,
                         const std::string       &units,
                         const TExpress::Options &o)
{
    o.writer->open("TransExpress_scatter.R");

    /*
    o.writer->write(RWriter::scatter(stats, ChrT, "????",
                                     "TransExpress",
                                     "Expected concentration (attomol/ul)",
                                     "Measured coverage (FPKM)",
                                     "Expected concentration (log2 attomol/ul)",
                                     "Measured coverage (log2 FPKM)"));
     */
    
    o.writer->close();
}

void TExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto m = std::map<TExpress::Metrics, std::string>
    {
        { TExpress::Metrics::Gene,    "gene"    },
        { TExpress::Metrics::Isoform, "isoform" },
    };
    
    const auto stats = TExpress::analyze(files, o);
    const auto units = m.at(o.metrs);
    
    /*
     * Generating summary statistics
     */

    o.info("Generating summary statistics");
    o.writer->open("TransExpress_summary.stats");
    
    for (auto i = 0; i < files.size(); i++)
    {
        o.writer->write(generateSummary(stats[i], files[i], units, o));

        if (i != files.size() - 1)
        {
            o.writer->write("\n\n");
        }
    }
    
    o.writer->close();
    
    /*
     * Generating summary statistics for each sample
     */
    
    std::vector<LinearStats> data;
    std::vector<MappingStats> data_;
    std::vector<SequinHist> hists;
    
    for (auto i = 0; i < files.size(); i++)
    {
        /*
         * Generating detailed statistics for the sequins
         */
        
        o.writer->open("TransExpress_quins.csv");
        o.writer->write(StatsWriter::writeCSV(stats[i].data, "Expected concentration (attomol/ul)", "Measured abundance (attomol/ul)"));
        o.writer->close();
        

        /*
         * Generating for AbundAbund
         */
        
        //o.writer->open("VarExpress_abundAbund.R");
        //o.writer->write(RWriter::createScript("VarExpress_quins.csv", PlotTAbundAbund()));
        //o.writer->close();

        // Generating scatter plot for the sample
        //writeScatter(stats[i], files[i], units, o);
        
        data_.push_back(stats[i]);
        hists.push_back(stats[i].hist);
        data.push_back(stats[i].data);
    }
    
    /*
     * Generating summary statistics for all samples
     */
    
    if (files.size() >= 2)
    {
        /*
         * Generating pooled summary statistics
         */
        
        o.writer->open("TransExpress_pooled.stats");
        o.writer->write(StatsWriter::inflectSummary(o.rChrT, o.rEndo, files, hists, data_, data, units));
        o.writer->close();
        
        /*
         * Generating pooled abundance
         */
        
        o.writer->open("TransExpress_pooled.R");
        o.writer->write(RWriter::scatterPool("TExpress_quins.csv"));
        o.writer->close();
    }

    /*
     * Generating spliced plot for all samples (but only if we have the isoforms...)
     */
    
    if (o.metrs == TExpress::Metrics::Isoform)
    {
        o.writer->open("TransExpress_splice.R");
        o.writer->write(RWriter::createSplice("TExpress_quins.csv"));
        o.writer->close();
    }
    
//    /*
//     * Generating CSV of expression for all samples
//     */
//    
//    const auto &r = Standard::instance().r_trans;
//    
//    o.writer->open("TExpress_quinsA.csv");
//    
//    switch (o.metrs)
//    {
//        case Metrics::Gene:
//        {
//            o.writer->write(writeSampleCSV(r.geneIDs(ChrT), stats, o));
//            break;
//        }
//            
//        case Metrics::Isoform:
//        {
//            o.writer->write(writeSampleCSV(r.seqIDs(), stats, o));
//            break;
//        }
//    }
//    
//    o.writer->close();
}