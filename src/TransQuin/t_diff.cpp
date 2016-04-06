/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include <ss/misc.hpp>
#include "TransQuin/t_diff.hpp"
#include "parsers/parser_edgeR.hpp"
#include "parsers/parser_sleuth.hpp"
#include "parsers/parser_DESeq2.hpp"
#include "parsers/parser_cdiffs.hpp"
#include "parsers/parser_HTSeqCount.hpp"

using namespace Anaquin;

typedef TDiffs::Metrics  Metrics;
typedef DiffTest::Status Status;
typedef TDiffs::Software Software;

std::vector<std::string> TDiffs::classify(const std::vector<double> &qs, const std::vector<double> &folds, double qCut, double foldCut)
{
    assert(qs.size() == folds.size());
    
    std::vector<std::string> r;
    
    for (auto i = 0; i < qs.size(); i++)
    {
        // Differential expressed
        if (qs[i] <= qCut)
        {
            r.push_back(fabs(folds[i]) <= foldCut ? "FP" : "TP");
        }
        
        // Non-differential expressed
        else
        {
            r.push_back(fabs(folds[i]) <= foldCut ? "TN" : "FN");
        }
    }
    
    return r;
}

template <typename T> void classifyChrT(TDiffs::Stats &stats, const T &t, const TDiffs::Options &o)
{
    assert(t.cID == ChrT);
    
    if (t.status != Status::Tested)
    {
        return;
    }
    
    const auto &id = t.id;
    const auto &r  = Standard::instance().r_trans;

    // Known fold change
    Fold known = NAN;
    
    // It's NAN if the sequin defined in reference but not in mixture
    Fold measured = NAN;

    switch (o.metrs)
    {
        case Metrics::Gene:
        {
            if (t.status == Status::Tested && stats.hist.count(t.id))
            {
                const auto *g = r.findGene(t.cID, id);
                
                if (g)
                {
                    stats.hist.at(id)++;

                    // Calculate the known fold-change between B and A
                    known = (g->abund(Mix_2) / g->abund(Mix_1));
                    
                    // This is not on the log scale, so it can't be zero...
                    assert(known);

                    // Measured fold-change between the two mixtures
                    measured = t.logF;
                    
                    // Turn it back to the original scale
                    measured = std::pow(2, measured);
                }

                stats.data[t.cID].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            }

            break;
        }

        case Metrics::Isoform:
        {
            if (t.status == Status::Tested && stats.hist.count(id))
            {
                const auto *seq = r.match(id);
                
                if (seq)
                {
                    // Known fold-change between the two mixtures
                    known = seq->abund(Mix_2) / seq->abund(Mix_1);
                    
                    stats.hist.at(id)++;
                    
                    // Measured fold-change between the two mixtures
                    measured = t.logF;
                    
                    // Turn it back to the original scale
                    measured = std::pow(2, measured);
                }
                
                stats.data[ChrT].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            }

            break;
        }
    }
    
    stats.data[t.cID].eLogFs.push_back(log2(known));
}

template <typename T> void classifyEndo(TDiffs::Stats &stats, const T &t, const TDiffs::Options &)
{
    stats.data[t.cID].eLogFs.push_back(NAN);
}

template <typename T> void update(TDiffs::Stats &stats, const T &t, const TDiffs::Options &o)
{
    if (t.cID == ChrT)
    {
        stats.n_chrT++;
        classifyChrT(stats, t, o);
    }
    else
    {
        stats.n_endo++;
        classifyEndo(stats, t, o);
    }

    stats.data[t.cID].ids.push_back(t.id);
    
    /*
     * The expected log-folds is done in the classifer (because it can only be done for synthetic)
     */
    
    if (t.status == Status::Tested)
    {
        stats.data[t.cID].ps.push_back(t.p);
        stats.data[t.cID].logFs.push_back(t.logF);
        stats.data[t.cID].logFSEs.push_back(t.logFSE);
        stats.data[t.cID].baseMeans.push_back(t.baseMean);
    }
    else
    {
        stats.data[t.cID].ps.push_back(NAN);
        stats.data[t.cID].logFs.push_back(NAN);
        stats.data[t.cID].logFSEs.push_back(NAN);
        stats.data[t.cID].baseMeans.push_back(NAN);
    }
}

template <typename Functor> TDiffs::Stats calculate(const TDiffs::Options &o, Functor f)
{
    const auto &r = Standard::instance().r_trans;

    TDiffs::Stats stats;

    stats.data[Endo];
    stats.data[ChrT];

    switch (o.metrs)
    {
        case Metrics::Gene:    { stats.hist = r.geneHist(ChrT); break; }
        case Metrics::Isoform: { stats.hist = r.hist();         break; }
    }

    assert(!stats.hist.empty());

    o.info("Parsing input files");
    f(stats);
    
    o.info("Calculating detection limit");

    switch (o.metrs)
    {
        case Metrics::Gene:    { stats.limit = r.absoluteGene(stats.hist); break; }
        //case Metrics::Isoform: { stats.limit = r.limitIsof(stats.hist); break; }
        default: { break; } // TODO: Please fix this
    }
    
    /*
     * We shouldn't assume any ordering in the inputs, we'll sort the data and assume
     * uniqueness in the genome.
     */
    
    // Used for quickly sort the count tables
    std::set<FeatureID> isChrT, isEndo;
    
    for (auto &i : stats.data)
    {
        const auto p = SS::sortPerm(i.second.ids, [&](const FeatureID &x, const FeatureID &y)
        {
            return x < y;
        });
        
        i.second.ps        = SS::applyPerm(i.second.ps,       p);
        i.second.ids       = SS::applyPerm(i.second.ids,      p);
        i.second.logFs     = SS::applyPerm(i.second.logFs,    p);
        i.second.eLogFs    = SS::applyPerm(i.second.eLogFs,   p);
        i.second.logFSEs   = SS::applyPerm(i.second.logFSEs,  p);
        i.second.baseMeans = SS::applyPerm(i.second.baseMeans, p);
        
        for (const auto &id : i.second.ids)
        {
            assert(!isChrT.count(id) && !isEndo.count(id));
            
            if (i.first == ChrT)
            {
                isChrT.insert(id);
            }
            else
            {
                isEndo.insert(id);
            }
        }
    }
    
    return stats;
}

TDiffs::Stats TDiffs::analyze(const std::vector<DiffTest> &tests, const Options &o)
{
    return calculate(o, [&](TDiffs::Stats &stats)
    {
        for (auto &test : tests)
        {
            update(stats, test, o);
        }
    });
}

TDiffs::Stats TDiffs::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](TDiffs::Stats &stats)
    {
        switch (o.dSoft)
        {
/*
            case Software::ParserSleuth:
            {
                ParserSleuth::parse(file, [&](const DiffTest &t, const ParserProgress &)
                {
                    update(stats, t, o);
                });

                break;
            }
*/
            case Software::DESeq2:
            {
                ParserDESeq2::parse(file, [&](const DiffTest &t, const ParserProgress &)
                {
                    update(stats, t, o);
                });

                break;                
            }
                
            case Software::edgeR:
            {
                break;
            }

            case Software::Cuffdiff:
            {
                ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
                {
                    update(stats, t, o);
                });

                break;
            }
        }
    });
}

static void writeDifferent(const FileName &file, const TDiffs::Stats &stats, const TDiffs::Options &o)
{
    /*
     * Generating a file for differential results. The file should list the relevant information for
     * plotting an MA and LODR plot.
     *
     * We will need the following information:
     *
     *     - BaseMean
     *     - Expected LF
     *     - LogFold
     *     - LogFold SE
     *     - PValue
     */
    
    o.writer->open(file);
    o.writer->write(",baseMean,elfc,lfc,lfcSE,pval");

    for (const auto &i : stats.data)
    {
        const auto &ps        = i.second.ps;
        const auto &ids       = i.second.ids;
        const auto &logFs     = i.second.logFs;
        const auto &eLogFs    = i.second.eLogFs;
        const auto &logFSEs   = i.second.logFSEs;
        const auto &baseMeans = i.second.baseMeans;

        for (auto j = 0; j < ids.size(); j++)
        {
            if (isnan(ps[j]))
            {
                o.writer->write((boost::format("%1%,NA,NA,NA") % ids[j]).str());
            }
            else
            {
                o.writer->write((boost::format("%1%,%2%,%3%,%4%,%5%,%6%") % ids[j]
                                                                          % toNA(baseMeans[j])
                                                                          % toNA(eLogFs[j])
                                                                          % toNA(logFs[j])
                                                                          % toNA(logFSEs[j])
                                                                          % toNA(ps[j])).str());
            }
        }
    }
    
    o.writer->close();
}

void TDiffs::report(const FileName &file, const Options &o)
{
    const auto stats = TDiffs::analyze(file, o);
    
    const auto m = std::map<Metrics, std::string>
    {
        { Metrics::Gene,    "gene"    },
        { Metrics::Isoform, "isoform" },
    };

    const auto units = m.at(o.metrs);
    
    o.info("Generating statistics");
    
    /*
     * There's no need to write summary for each replicate because differential analysis has already incorporated them.
     */
    
    /*
     *  - Summary statistics
     *  - Sequin statistics
     *  - Scatter plot
     *  - ROC plot
     *  - MA plot
     *  - LODR Plot
     */

    /*
     * Generating summary statistics
     */
    
    o.writer->open("TransDiff_summary.stats");
    //o.writer->write(StatsWriter::linear(file, stats, ChrT, units));
    o.writer->close();
    
    /*
     * Generating scatter plot for the log-fold changes
     */
    
    o.writer->open("TransDiff_scatter.R");
    o.writer->write(RWriter::scatter(stats, ChrT, "????", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
    o.writer->close();

    /*
     * Generating ROC plot
     */
    
    o.writer->open("TransDiff_ROC.R");
    //o.writer->write(RWriter::createROC_T(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps, units));
    o.writer->close();

    /*
     * Generating differential results (CSV)
     */
    
    writeDifferent("TransDiff_diffs.csv", stats, o);

    /*
     * Generating MA plot
     */
    
    o.writer->open("TransDiff_MA.R");
    o.writer->write(RWriter::createMA("TransDiff_diffs.csv", units));
    o.writer->close();
    
    /*
     * Generating LODR plot
     */
    
    o.writer->open("TransDiff_LODR.R");
    o.writer->write(RWriter::createLODR_T("TransDiff_diffs.csv"));
    o.writer->close();
}
