/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "RnaQuin/r_fold.hpp"
#include "parsers/parser_diff.hpp"

using namespace Anaquin;

extern Scripts PlotTROC();

typedef RFold::Metrics Metrics;

std::vector<std::string> RFold::classify(const std::vector<double> &qs, const std::vector<double> &folds, double qCut, double foldCut)
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

template <typename T> void classifySyn(RFold::Stats &stats, const T &t, const RFold::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    auto f = [&](const SequinID &id, Concent exp)
    {
        stats.data[id].p    = t.p;
        stats.data[id].q    = t.q;
        stats.data[id].exp  = exp;
        stats.data[id].obs  = t.logF;
        stats.data[id].se   = t.logFSE;
        stats.data[id].mean = t.mean;
        
        if (!isnan(exp) && !isnan(t.logF))
        {
            if (isnan(stats.limit.abund) || exp < stats.limit.abund)
            {
                stats.limit.id = id;
                stats.limit.abund = exp;
            }

            stats.add(id, exp, t.logF);
        }
    };
    
    switch (o.metrs)
    {
        case Metrics::Gene:
        {
            const auto match = r.findGene(t.cID, t.id);
            
            if (match)
            {
                stats.hist.at(t.cID).at(t.id)++;
                
                const auto exp_1 = r.concent(t.id, Mix_1);
                const auto exp_2 = r.concent(t.id, Mix_2);
                
                f(t.id, log2(exp_2 / exp_1));
            }
            else
            {
                o.warn(t.id + " not found");
            }

            break;
        }

        case Metrics::Isoform:
        {
            const auto match = r.match(t.id);
            
            if (match)
            {
                stats.hist.at(t.cID).at(t.id)++;
                
                const auto e1 = match->concent(Mix_1);
                const auto e2 = match->concent(Mix_2);
                
                f(t.id, log2(e2 / e1));
            }
            else
            {
                o.warn(t.id + " not found");
            }
            
            break;
        }
    }
}

template <typename T> void update(RFold::Stats &stats, const T &x, const RFold::Options &o)
{
    typedef DiffTest::Status Status;
    
    if (Standard::isSynthetic(x.cID))
    {
        stats.n_syn++;
        classifySyn(stats, x, o);
    }
    else
    {
        stats.n_gen++;
    }
}

template <typename Functor> RFold::Stats calculate(const RFold::Options &o, Functor f)
{
    const auto &r = Standard::instance().r_trans;

    RFold::Stats stats;

    switch (o.metrs)
    {
        case Metrics::Gene:    { stats.hist = r.histGene(); break; }
        case Metrics::Isoform: { stats.hist = r.histIsof(); break; }
    }

    assert(!stats.hist.empty());

    o.info("Parsing input files");
    f(stats);
    
    return stats;
}

RFold::Stats RFold::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](RFold::Stats &stats)
    {
        ParserDiff::parse(file, [&](const ParserDiff::Data &x, const ParserProgress &)
        {
            update(stats, x, o);
        });
    });
}

void RFold::report(const FileName &file, const Options &o)
{
    const auto m = std::map<Metrics, std::string>
    {
        { Metrics::Gene,    "genes"    },
        { Metrics::Isoform, "isoforms" },
    };

    switch (o.metrs)
    {
        case Metrics::Gene:    { o.info("Gene Differential");    break; }
        case Metrics::Isoform: { o.info("Isoform Differential"); break; }
    }
    
    const auto stats = RFold::analyze(file, o);
    const auto units = m.at(o.metrs);
    
    o.info("Generating statistics");
    
    // Eg: DESeq2
    const auto shouldLOD = true;
    
    /*
     * Generating RnaFoldChange_summary.stats
     */

    RFold::generateSummary("RnaFoldChange_summary.stats", file, stats, o, units);

    /*
     * Generating RnaFoldChange_sequins.csv
     */

    RFold::generateCSV("RnaFoldChange_sequins.csv", stats, o);
    
    /*
     * Generating RnaFoldChange_fold.R
     */
    
    o.generate("RnaFoldChange_fold.R");
    o.writer->open("RnaFoldChange_fold.R");
    o.writer->write(RWriter::createFold("RnaFoldChange_sequins.csv",
                                        "Fold Change",
                                        "Expected fold change (log2)",
                                        "Measured fold change (log2)",
                                        "Expected",
                                        "Measured", false));
    o.writer->close();

    /*
     * Generating RnaFoldChange_ROC.R
     */
    
    o.generate("RnaFoldChange_ROC.R");
    o.writer->open("RnaFoldChange_ROC.R");
    o.writer->write(RWriter::createScript("RnaFoldChange_sequins.csv", PlotTROC()));
    o.writer->close();

    /*
     * Generating RnaFoldChange_LODR.R
     */
    
    if (shouldLOD)
    {
        RFold::generateLODR("RnaFoldChange_LODR.R", "RnaFoldChange_sequins.csv", o);
    }
}