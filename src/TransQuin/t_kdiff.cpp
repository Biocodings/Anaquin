/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include <ss/misc.hpp>
#include "TransQuin/t_kdiff.hpp"
#include "parsers/parser_sleuth.hpp"

using namespace Anaquin;

typedef TKDiff::Metrics  Metrics;
typedef DiffTest::Status Status;
typedef TKDiff::Software Software;

template <typename T> void update(TKDiff::Stats &stats, const T &t, const TKDiff::Options &o)
{
//    if (t.cID == ChrT)
//    {
//        stats.n_chrT++;
//        classifyChrT(stats, t, o);
//    }
//    else
//    {
//        stats.n_endo++;
//        classifyEndo(stats, t, o);
//    }
//
//    stats.data[t.cID].ids.push_back(t.id);
//    
//    /*
//     * The expected log-folds is done in the classifer (because it can only be done for synthetic)
//     */
//    
//    if (t.status == Status::Tested)
//    {
//        stats.data[t.cID].ps.push_back(t.p);
//        stats.data[t.cID].logFs.push_back(t.logF);
//        stats.data[t.cID].logFSEs.push_back(t.logFSE);
//        stats.data[t.cID].baseMeans.push_back(t.baseMean);
//    }
//    else
//    {
//        stats.data[t.cID].ps.push_back(NAN);
//        stats.data[t.cID].logFs.push_back(NAN);
//        stats.data[t.cID].logFSEs.push_back(NAN);
//        stats.data[t.cID].baseMeans.push_back(NAN);
//    }
}

template <typename Functor> TKDiff::Stats calculate(const TKDiff::Options &o, Functor f)
{
    //const auto &r = Standard::instance().r_trans;

    TKDiff::Stats stats;

//    //stats.data[Endo];
//    //stats.data[ChrT];
//
//    switch (o.metrs)
//    {
//        case Metrics::Gene:    { stats.hist = r.geneHist(ChrT); break; }
//        case Metrics::Isoform: { stats.hist = r.hist();         break; }
//    }
//
//    assert(!stats.hist.empty());
//
//    o.info("Parsing input files");
//    f(stats);
//    
//    o.info("Calculating detection limit");
//
//    switch (o.metrs)
//    {
//        case Metrics::Gene:    { stats.limit = r.limitGene(stats.hist); break; }
//        //case Metrics::Isoform: { stats.limit = r.limitIsof(stats.hist); break; }
//        default: { break; } // TODO: Please fix this
//    }
//    
//    /*
//     * We shouldn't assume any ordering in the inputs, we'll sort the data and assume
//     * uniqueness in the genome.
//     */
//    
//    // Used for quickly sort the count tables
//    std::set<FeatureID> isChrT, isEndo;
//    
//    for (auto &i : stats.data)
//    {
//        const auto p = SS::sortPerm(i.second.ids, [&](const FeatureID &x, const FeatureID &y)
//        {
//            return x < y;
//        });
//        
//        i.second.ps        = SS::applyPerm(i.second.ps,       p);
//        i.second.ids       = SS::applyPerm(i.second.ids,      p);
//        i.second.logFs     = SS::applyPerm(i.second.logFs,    p);
//        i.second.eLogFs    = SS::applyPerm(i.second.eLogFs,   p);
//        i.second.logFSEs   = SS::applyPerm(i.second.logFSEs,  p);
//        i.second.baseMeans = SS::applyPerm(i.second.baseMeans, p);
//        
//        for (const auto &id : i.second.ids)
//        {
//            assert(!isChrT.count(id) && !isEndo.count(id));
//            
//            if (i.first == ChrT)
//            {
//                isChrT.insert(id);
//            }
//            else
//            {
//                isEndo.insert(id);
//            }
//        }
//    }
//    
//    if (stats.counts)
//    {
//        readCounts(isChrT, isEndo, stats, o);
//    }
    
    return stats;
}

TKDiff::Stats TKDiff::analyze(const FileName &file, const Options &o)
{
//    return calculate(o, [&](TKDiff::Stats &stats)
//    {
//        switch (o.dSoft)
//        {
//            case DiffSoft::ParserSleuth:
//            {
//                ParserSleuth::parse(file, [&](const DiffTest &t, const ParserProgress &)
//                {
//                    update(stats, t, o);
//                });
//
//                break;
//            }
//
//            case DiffSoft::DESeq2:
//            {
//                ParserDESeq2::parse(file, [&](const DiffTest &t, const ParserProgress &)
//                {
//                    update(stats, t, o);
//                });
//
//                break;                
//            }
//                
//            case DiffSoft::edgeR:
//            {
//                break;
//            }
//
//            case DiffSoft::Cuffdiff:
//            {
//                ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
//                {
//                    update(stats, t, o);
//                });
//
//                break;
//            }
//        }
//    });
    
    throw "";
}

//static void writeDifferent(const FileName &file, const TKDiff::Stats &stats, const TKDiff::Options &o)
//{
//    /*
//     * Generating a file for differential results. The file should list the relevant information for
//     * plotting an MA and LODR plot.
//     *
//     * We will need the following information:
//     *
//     *     - BaseMean
//     *     - Expected LF
//     *     - LogFold
//     *     - LogFold SE
//     *     - PValue
//     */
//    
//    o.writer->open(file);
//    o.writer->write(",baseMean,elfc,lfc,lfcSE,pval");
//
//    for (const auto &i : stats.data)
//    {
//        const auto &ps        = i.second.ps;
//        const auto &ids       = i.second.ids;
//        const auto &logFs     = i.second.logFs;
//        const auto &eLogFs    = i.second.eLogFs;
//        const auto &logFSEs   = i.second.logFSEs;
//        const auto &baseMeans = i.second.baseMeans;
//
//        for (auto j = 0; j < ids.size(); j++)
//        {
//            if (isnan(ps[j]))
//            {
//                o.writer->write((boost::format("%1%,NA,NA,NA") % ids[j]).str());
//            }
//            else
//            {
//                o.writer->write((boost::format("%1%,%2%,%3%,%4%,%5%,%6%") % ids[j]
//                                                                          % toNA(baseMeans[j])
//                                                                          % toNA(eLogFs[j])
//                                                                          % toNA(logFs[j])
//                                                                          % toNA(logFSEs[j])
//                                                                          % toNA(ps[j])).str());
//            }
//        }
//    }
//    
//    o.writer->close();
//}

void TKDiff::report(const FileName &file, const Options &o)
{
//    const auto stats = TKDiff::analyze(file, o);
//    
//    const auto m = std::map<Metrics, std::string>
//    {
//        { Metrics::Gene,    "gene"    },
//        { Metrics::Isoform, "isoform" },
//    };
//
//    const auto units = m.at(o.metrs);
//    
//    o.info("Generating statistics");
//    
//    /*
//     * There's no need to write summary for each replicate because differential analysis has already incorporated them.
//     */
//    
//    /*
//     *  - Summary statistics
//     *  - Sequin statistics
//     *  - Scatter plot
//     *  - ROC plot
//     *  - MA plot
//     *  - LODR Plot
//     */
//
//    /*
//     * Generating summary statistics
//     */
//    
//    o.writer->open("TransDiffs_summary.stats");
//    //o.writer->write(StatsWriter::linear(file, stats, ChrT, units));
//    o.writer->close();
//    
//    /*
//     * Generating scatter plot for the log-fold changes
//     */
//    
//    o.writer->open("TransDiffs_scatter.R");
//    o.writer->write(RWriter::scatter(stats, ChrT, "????", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
//    o.writer->close();
//
//    /*
//     * Generating ROC plot
//     */
//    
//    o.writer->open("TransDiffs_ROC.R");
//    //o.writer->write(RWriter::createROC_T(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps, units));
//    o.writer->close();
//
//    /*
//     * Generating differential results (CSV)
//     */
//    
//    writeDifferent("TransDiffs_diffs.csv", stats, o);
//
//    /*
//     * Generating MA plot
//     */
//    
//    o.writer->open("TransDiffs_MA.R");
//    o.writer->write(RWriter::createMA("TransDiffs_diffs.csv", units));
//    o.writer->close();
//    
//    /*
//     * Generating LODR plot
//     */
//    
//    o.writer->open("TransDiffs_LODR.R");
//    o.writer->write(RWriter::createLODR_T("TransDiffs_diffs.csv"));
//    o.writer->close();
}
