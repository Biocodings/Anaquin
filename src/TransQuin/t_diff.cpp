/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "TransQuin/t_diff.hpp"
#include "parsers/parser_edgeR.hpp"
#include "parsers/parser_sleuth.hpp"
#include "parsers/parser_DESeq2.hpp"
#include "parsers/parser_cdiff.hpp"
#include "parsers/parser_HTSeqCount.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotFold();

// Defined in resources.cpp
extern Scripts PlotLODR();

// Defined in resources.cpp
extern Scripts PlotTROC();

// Defined in resources.cpp
extern Scripts PlotMA();

typedef TDiff::Metrics  Metrics;
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
            if (stats.hist.count(t.id))
            {
                const auto *match = r.findGene(t.cID, id);
                
                if (match)
                {
                    stats.hist.at(id)++;

                    // Calculate the known fold-change between B and A
                    known = match->abund(Mix_2) / match->abund(Mix_1);
                    
                    // This is not on the log scale, so it can't be zero...
                    assert(known);

                    // Measured fold-change between the two mixtures
                    measured = t.logF;
                    
                    // Turn it back to the original scale
                    measured = std::pow(2, measured);
                }

                stats.add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            }

            break;
        }

        case Metrics::Isoform:
        {
            if (stats.hist.count(id))
            {
                const auto *match = r.match(id);
                
                if (match)
                {
                    stats.hist.at(id)++;

                    // Known fold-change between the two mixtures
                    known = match->abund(Mix_2) / match->abund(Mix_1);
                    
                    // Measured fold-change between the two mixtures
                    measured = t.logF;
                    
                    // Turn it back to the original scale
                    measured = std::pow(2, measured);
                }
                
                stats.add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            }

            break;
        }
    }
    
    stats.eLogFs.push_back(log2(known));
}

template <typename T> void update(TDiff::Stats &stats, const T &t, const TDiff::Options &o)
{
    typedef DiffTest::Status Status;
    
    if (t.status != Status::Tested)
    {
        return;
    }
    else if (t.cID == ChrT)
    {
        stats.n_chrT++;
        classifyChrT(stats, t, o);
    }
    else
    {
        stats.n_endo++;
        stats.eLogFs.push_back(NAN);
    }

    /*
     * The expected log-fold is done in the classifer because it can only happen for synthetic
     */

    stats.ids.push_back(t.id);
    
    if (t.status == Status::Tested)
    {
        stats.ps.push_back(t.p);
        stats.qs.push_back(t.q);
        stats.mLogFs.push_back(t.logF);
        stats.logFSEs.push_back(t.logFSE);
        stats.baseMeans.push_back(t.baseMean);
    }
    else
    {
        stats.ps.push_back(NAN);
        stats.qs.push_back(NAN);
        stats.mLogFs.push_back(NAN);
        stats.logFSEs.push_back(NAN);
        stats.baseMeans.push_back(NAN);
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
        case Metrics::Gene:    { stats.limit = r.absoluteGene(stats.hist); break; }
        //case Metrics::Isoform: { stats.limit = r.limitIsof(stats.hist); break; }
        default: { break; }
    }

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
            case Software::Sleuth:
            {
                ParserSleuth::parse(file, [&](const ParserSleuth::Data &data, const ParserProgress &)
                {
                    update(stats, data, o);
                });

                break;
            }
                
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
                ParserEdgeR::parse(file, [&](const DiffTest &t, const ParserProgress &)
                {
                    update(stats, t, o);
                });

                break;
            }

            case Software::Cuffdiff:
            {
                ParserCDiff::parse(file, [&](const ParserCDiff::Data &data, const ParserProgress &)
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
    o.writer->write(StatsWriter::linearSummary(file, o.rChrT, stats, stats.hist));
    o.writer->close();
    
    /*
     * 2. Generating differential results
     */
    
    o.info("Generating TransDiff_quins.csv");
    o.writer->open("TransDiff_quins.csv");
    o.writer->write(TDiff::writeCSV(stats, o));
    o.writer->close();
    
    /*
     * 3. Generating log-fold plot
     */
    
    o.info("Generating TransDiff_fold.R");
    o.writer->open("TransDiff_fold.R");
    o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotFold()));
    o.writer->close();

    /*
     * 4. Generating ROC plot
     */
    
    o.info("Generating TransDiff_ROC.R");
    o.writer->open("TransDiff_ROC.R");
    o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotTROC()));
    o.writer->close();

    /*
     * 5. Generating LODR plot
     */
    
    o.info("Generating TransDiff_LODR.R");
    o.writer->open("TransDiff_LODR.R");
    o.writer->write(RWriter::createScript("TransDiff_quins.csv", PlotLODR()));
    o.writer->close();
    
    /*
     * 6. Generating MA plot
     */
    
    if (!o.counts.empty())
    {
        o.info("Generating TransDiff_MA.R");
        o.writer->open("TransDiff_MA.R");
        o.writer->write(RWriter::createScript("TransDiff_counts.csv", PlotMA()));
        o.writer->close();
    }
}
