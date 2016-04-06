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

// Defined in resources.cpp
extern Scripts PlotFold();

// Defined in resources.cpp
extern Scripts PlotROC_T();

// Defined in resources.cpp
extern Scripts PlotMA();

// Defined in resources.cpp
extern Scripts PlotLODR_T();

typedef TDiff::Metrics  Metrics;
typedef DiffTest::Status Status;
typedef TDiff::Software Software;

std::vector<std::string> TDiff::classify(const std::vector<double> &qs, const std::vector<double> &folds, double qCut, double foldCut)
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

template <typename T> void classifyChrT(TDiff::Stats &stats, const T &t, const TDiff::Options &o)
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

                stats.data.add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
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
                
                stats.data.add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            }

            break;
        }
    }
    
    stats.data.eLogFs.push_back(log2(known));
}

template <typename T> void classifyEndo(TDiff::Stats &stats, const T &t, const TDiff::Options &)
{
    stats.data.eLogFs.push_back(NAN);
}

template <typename T> void update(TDiff::Stats &stats, const T &t, const TDiff::Options &o)
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

    stats.data.ids.push_back(t.id);
    
    /*
     * The expected log-folds is done in the classifer (because it can only be done for synthetic)
     */
    
    if (t.status == Status::Tested)
    {
        stats.data.ps.push_back(t.p);
        stats.data.logFs.push_back(t.logF);
        stats.data.logFSEs.push_back(t.logFSE);
        stats.data.baseMeans.push_back(t.baseMean);
    }
    else
    {
        stats.data.ps.push_back(NAN);
        stats.data.logFs.push_back(NAN);
        stats.data.logFSEs.push_back(NAN);
        stats.data.baseMeans.push_back(NAN);
    }
}

template <typename Functor> TDiff::Stats calculate(const TDiff::Options &o, Functor f)
{
    const auto &r = Standard::instance().r_trans;

    TDiff::Stats stats;

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
        case Metrics::Gene:    { stats.data.limit = r.absoluteGene(stats.hist); break; }
        //case Metrics::Isoform: { stats.limit = r.limitIsof(stats.hist); break; }
        default: { break; } // TODO: Please fix this
    }
    
    /*
     * We shouldn't assume any ordering in the inputs, we'll sort the data and assume
     * uniqueness in the genome.
     */
    
    // Used for quickly sort the count tables
//    std::set<FeatureID> isChrT, isEndo;
//    
//    for (auto &i : stats)
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
    
    return stats;
}

TDiff::Stats TDiff::analyze(const std::vector<DiffTest> &tests, const Options &o)
{
    return calculate(o, [&](TDiff::Stats &stats)
    {
        for (auto &test : tests)
        {
            update(stats, test, o);
        }
    });
}

TDiff::Stats TDiff::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](TDiff::Stats &stats)
    {
        switch (o.dSoft)
        {
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
                ParserCDiffs::parse(file, [&](const ParserCDiffs::Data &data, const ParserProgress &)
                {
                    update(stats, data, o);
                });

                break;
            }
        }
    });
}

void TDiff::report(const FileName &file, const Options &o)
{
    const auto m = std::map<Metrics, std::string>
    {
        { Metrics::Gene,    "gene"    },
        { Metrics::Isoform, "isoform" },
    };

    const auto stats = TDiff::analyze(file, o);
    const auto units = m.at(o.metrs);
    
    o.info("Generating statistics");
    
    /*
     * 1. Generating summary statistics
     */
    
    o.info("Generating TransDiff_summary.stats");
    o.writer->open("TransDiff_summary.stats");
    o.writer->write(StatsWriter::linearSummary(file, "", stats.data, stats.hist));
    o.writer->close();
    
    /*
     * 2. Generating differential results (CSV)
     */
    
    o.info("Generating TransDiff_quins.csv");
    o.writer->open("TransDiff_quins.csv");
    o.writer->write(TDiff::writeCSV(stats, o));
    o.writer->close();
    
    /*
     * 3. Generating scatter plot for the log-fold changes
     */
    
    o.info("Generating TransDiff_fold.R");
    o.writer->open("TransDiff_fold.R");
    o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotFold()));
    //o.writer->write(RWriter::scatter(stats, ChrT, "????", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
    o.writer->close();

    /*
     * 4. Generating ROC plot
     */
    
    o.info("Generating TransDiff_ROC.R");
    o.writer->open("TransDiff_ROC.R");
    //o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotROC_T()));
    //o.writer->write(RWriter::createROC_T(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps, units));
    o.writer->close();

    /*
     * 5. Generating MA plot
     */
    
    if (!o.counts.empty())
    {
        o.info("Generating TransDiff_MA.R");
        o.writer->open("TransDiff_MA.R");
        o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotMA()));
        //o.writer->write(RWriter::createMA("TransDiff_diffs.csv", units));
        o.writer->close();
    }
    
    /*
     * 6. Generating LODR plot
     */
    
    o.info("Generating TransDiff_LODR.R");
    o.writer->open("TransDiff_LODR.R");
    //o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotLODR_T()));
    //o.writer->write(RWriter::createLODR_T("TransDiff_diffs.csv"));
    o.writer->close();
}
