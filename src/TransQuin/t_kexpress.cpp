#include <stdexcept>
#include "writers/r_writer.hpp"
#include "TransQuin/t_kexpress.hpp"
#include "parsers/parser_kallisto.hpp"
#include "parsers/parser_stringtie.hpp"
#include <ss/regression/segmented.hpp>

using namespace Anaquin;

typedef TKExpress::Metrics  Metrics;
typedef TKExpress::Software Software;

struct InternalKallistoData : public ParserKallisto::Data
{
    ChrID cID = ChrT;
    
    // Dummy value...
    Locus l;
};

template <typename T> void update(TKExpress::Stats &stats, const T &t, const TKExpress::Options &o)
{
    //    if (t.cID != ChrT)
    //    {
    //        stats.n_endo++;
    //    }
    //    else
    //    {
    //        stats.n_chrT++;
    //    }
    //
    //    if (t.cID == ChrT)
    //    {
    //        const auto &r = Standard::instance().r_trans;
    //
    //        switch (o.metrs)
    //        {
    //            case Metrics::Isoform:
    //            {
    //                const TransData *m = nullptr;
    //
    //                // Try to match by name if possible
    //                m = r.match(t.id);
    //
    //                if (!m)
    //                {
    //                    // Try to match by locus (de-novo assembly)
    //                    m = r.match(t.l, Overlap);
    //                }
    //
    //                if (!m)
    //                {
    //                    o.logWarn((boost::format("%1% not found. Unknown isoform.") % t.id).str());
    //                }
    //                else
    //                {
    //                    stats.hist.at(m->id)++;
    //
    //                    if (t.abund)
    //                    {
    //                        stats.data[t.cID].add(t.id, m->abund(Mix_1), t.abund);
    //                    }
    //                }
    //
    //                break;
    //            }
    //
    //            case Metrics::Gene:
    //            {
    //                const TransRef::GeneData *m = nullptr;
    //
    //                // Try to match by name if possible
    //                m = r.findGene(t.cID, t.id);
    //
    //                if (!m)
    //                {
    //                    // Try to match by locus (de-novo assembly)
    //                    m = r.findGene(t.cID, t.l, Contains);
    //                }
    //
    //                if (m)
    //                {
    //                    stats.hist.at(m->id)++;
    //
    //                    if (t.abund)
    //                    {
    //                        stats.data[t.cID].add(t.id, m->abund(Mix_1), t.abund);
    //                    }
    //                }
    //                else
    //                {
    //                    o.logWarn((boost::format("%1% not found. Unknown gene.") % t.id).str());
    //                }
    //
    //                break;
    //            }
    //        }
    //    }
}

TKExpress::Stats TKExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);
    
    //    return calculate(o.exp->fileToSample(file), o, [&](TKExpress::Stats &stats)
    //    {
    //        switch (o.soft)
    //        {
    //            case Software::Kallisto:
    //            {
    //                ParserKallisto::parse(file, [&](const ParserKallisto::Data &data, const ParserProgress &p)
    //                {
    //                    InternalKallistoData x;
    //
    //                    x.id    = data.id;
    //                    x.abund = data.abund;
    //
    //                    update(stats, x, o);
    //                });
    //
    //                break;
    //            }
    //        }
    //    });
    
    throw "";
}

//static void writeSummary(const TKExpress::Stats &stats,
//                         const FileName        &file,
//                         const std::string     &name,
//                         const Units           &units,
//                         const TKExpress::Options &o)
//{
////    o.writer->create(name);
////    o.writer->open(name + "/TransExpress_summary.stats");
////    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
////                                                o.rEndo,
////                                                std::vector<FileName>     { file  },
////                                                std::vector<SequinHist>   { stats.hist },
////                                                std::vector<MappingStats> { stats },
////                                                std::vector<LinearStats>  { stats.data.at(ChrT) },
////                                                units));
////    o.writer->close();
//}
//
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

void TKExpress::report(const std::vector<FileName> &files, const Options &o)
{
    //    const auto stats = TKExpress::analyze(files, o);
    //
    //    const auto m = std::map<TKExpress::Metrics, std::string>
    //    {
    //        { TKExpress::Metrics::Gene,    "gene"    },
    //        { TKExpress::Metrics::Isoform, "isoform" },
    //    };
    //
    //    const auto units = m.at(o.metrs);
    //
    //    o.info("Generating statistics");
    //
    //    /*
    //     * Generating summary statistics for each sample
    //     */
    //
    //    std::vector<TKExpress::Stats::Data> data;
    //    std::vector<MappingStats> data_;
    //    std::vector<SequinHist> hists;
    //
    //    for (auto i = 0; i < files.size(); i++)
    //    {
    //        // Generating summary statistics for the sample
    //        writeSummary(stats[i], files[i], o.exp->names().at(i), units, o);
    //
    //        o.writer->open(o.exp->names().at(i) + "/TransExpress_quins.csv");
    //        o.writer->write(StatsWriter::writeCSV(stats[i].data.at(ChrT), "Expected concentration (attomol/ul)", "Measured abundance (attomol/ul)"));
    //        o.writer->close();
    //
    //        // Generating scatter plot for the sample
    //        writeScatter(stats[i], files[i], o.exp->names().at(i), units, o);
    //
    //        data_.push_back(stats[i]);
    //        hists.push_back(stats[i].hist);
    //        data.push_back(stats[i].data.at(ChrT));
    //    }
    //
    //    /*
    //     * Generating CSV of expression for all samples
    //     */
    //
    ////    const auto &r = Standard::instance().r_trans;
    ////
    ////    o.writer->open("ABCD");
    ////
    ////    switch (o.lvl)
    ////    {
    ////        case Level::Gene:
    ////        {
    ////            o.writer->write(writeSampleCSV(r.geneIDs(ChrT), stats, o));
    ////            break;
    ////        }
    ////
    ////        case Level::Isoform:
    ////        {
    ////            o.writer->write(writeSampleCSV(r.seqIDs(), stats, o));
    ////            break;
    ////        }
    ////
    ////        case Level::Exon:
    ////        {
    ////            break;
    ////        }
    ////    }
    ////
    ////    o.writer->close();
    //
    //    /*
    //     * Generating a CSV of expression for all samples
    //     */
    //
    //    writeFPKM("TKExpress_FPKM.csv", stats, o);
    //
    //    /*
    //     * Generating summary statistics for all samples
    //     */
    //
    //    o.writer->open("TransExpress_pooled.stats");
    //    o.writer->write(StatsWriter::inflectSummary(o.rChrT, o.rEndo, files, hists, data_, data, units));
    //    o.writer->close();
    //
    //    /*
    //     * Generating scatter plot for all samples
    //     */
    //    
    //    o.writer->open("TransExpress_pooled.R");
    //    o.writer->write(RWriter::scatterPool("TKExpress_FPKM.csv"));
    //    o.writer->close();
    //
    //    /*
    //     * Generating spliced plot for all samples (but only if we have the isoforms...)
    //     */
    //
    //    if (o.metrs == TKExpress::Metrics::Isoform)
    //    {
    //        o.writer->open("TransExpress_Splice.R");
    //        o.writer->write(RWriter::createSplice("TKExpress_FPKM.csv"));
    //        o.writer->close();
    //    }
}