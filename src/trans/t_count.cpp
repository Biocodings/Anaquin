/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "trans/t_count.hpp"
#include "parsers/parser_cdiffs.hpp"

using namespace Anaquin;

typedef DiffTest::Status Status;
typedef TCount::Metrics  Metrics;
typedef TCount::Software Software;

std::vector<std::string> TCount::classify(const std::vector<double> &qs, const std::vector<double> &folds, double qCut, double foldCut)
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

template <typename T> void classifyChrT(TCount::Stats &stats, const T &t, const GenericID &id, Metrics metrs)
{
    assert(t.cID == ChrT);
    
    const auto &r = Standard::instance().r_trans;
    
    // Known fold change
    Fold known = NAN;
    
    // It's NAN if the sequin defined in reference but not in mixture
    Fold measured = NAN;
    
    /*
     * Differential expression at the gene level
     */
    
    auto g = [&](const GeneID &id, double fpkm_1, double fpkm_2)
    {
        const auto *g = r.findGene(t.cID, id);
        
        if (g)
        {
            // Calculate the known fold-change between B and A
            known = (g->abund(Mix_2) / g->abund(Mix_1));
        }
        
        if (g && !isnan(fpkm_1) && !isnan(fpkm_2) && fpkm_1 && fpkm_2)
        {
            stats.hist.at(id)++;
            
            // Measured fold-change between the two mixtures
            measured = fpkm_2 / fpkm_1;
        }
        
        stats.data[t.cID].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
    };
    
    switch (metrs)
    {
        case Metrics::Gene:
        {
            if ((t.status != Status::NotTested) && stats.hist.count(t.id))
            {
                g(t.id, t.fpkm_1, t.fpkm_2);
            }
            
            break;
        }
            
        case Metrics::Isoform:
        {
            if ((t.status == Status::NotTested) || !stats.hist.count(id))
            {
                return;
            }
            
            const auto *seq = r.match(id);
            
            if (seq)
            {
                // Known fold-change between the two mixtures
                known = seq->abund(Mix_2) / seq->abund(Mix_1);
            }
            
            if ((t.status != Status::NotTested) && t.fpkm_1 && t.fpkm_2)
            {
                stats.hist.at(id)++;
                
                // Measured fold-change between the two mixtures
                measured = t.fpkm_2 / t.fpkm_1;
            }
            
            stats.data[ChrT].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            
            break;
        }
    }

    stats.data[ChrT].seqs.push_back(id);
    stats.data[ChrT].ps.push_back(t.p);
    stats.data[ChrT].qs.push_back(t.q);
    stats.data[ChrT].logFCs.push_back(log2(known));
}

template <typename T> void classifyEndT(TCount::Stats &, const T &, const GenericID &, Metrics)
{
    /*
     * Obviously, we can't compare fold-changes for endogenous data. There's nothing else to do here...
     */
}

template <typename T> void update(TCount::Stats &stats, const T &t, const GenericID &id, Metrics metrs)
{
    if (t.cID == ChrT)
    {
        stats.n_chrT++;
        classifyChrT(stats, t, id, metrs);
    }
    else
    {
        stats.n_endo++;
        classifyEndT(stats, t, id, metrs);
    }
}

template <typename Functor> TCount::Stats calculate(const TCount::Options &o, Functor f)
{
    TCount::Stats stats;

    const auto &r = Standard::instance().r_trans;
    const auto cIDs = r.chromoIDs();

    /*
     * Initalize data for the synthetic chromosome (not needed for anything else).
     */

    stats.data[ChrT];

    const auto isoform = (o.metrs == Metrics::Isoform);
    o.logInfo(isoform ? "Isoform metrics" : "Gene metrics");
    
    // Construct for a histogram for the appropriate metrics
    stats.hist = isoform ? r.hist() : r.geneHist(ChrT);
    
    o.info("Parsing input file");

    f(stats);
    
    /*
     * Calculating detection limit
     */
    
    o.info("Calculating detection limit");
    
    stats.limit = isoform ? r.limit(stats.hist) : r.limitGene(stats.hist);
    
    return stats;
}

TCount::Stats TCount::analyze(const std::vector<DiffTest> &tests, const Options &o)
{
    return calculate(o, [&](TCount::Stats &stats)
    {
        for (auto &test : tests)
        {
            update(stats, test, test.id, o.metrs);
        }
    });
}

TCount::Stats TCount::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](TCount::Stats &stats)
    {
        switch (o.soft)
        {
            case Software::Cuffdiffs:
            {
                ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
                {
                    update(stats, t, o.metrs == Metrics::Isoform ? t.testID : t.id, o.metrs);
                });

                break;
            }
        }
    });
}

void TCount::report(const FileName &file, const Options &o)
{
    const auto stats = TCount::analyze(file, o);
    const auto units = (o.metrs == Metrics::Isoform) ? "isoforms" : "genes";
    
    o.info("Generating statistics");
    
    /*
     * There's no need to write summary for each replicate because differential analysis has already incorporated them.
     */
    
    /*
     * Synthetic
     * ---------
     *
     *    - Summary statistics
     *    - Sequin statistics
     *    - Scatter plot
     *    - ROC plot
     *    - MA plot
     *    - LODR Plot
     */

    /*
     * Endogenous
     * ----------
     *
     *    - MA Plot (merged with synthetic)
     */
    
    /*
     * Generating summary statistics
     */
    
    o.writer->open("TransDiffs_summary.stats");
    o.writer->write(StatsWriter::linear(file, stats, ChrT, units));
    o.writer->close();
    
    /*
     * Generating scatter plot
     */
    
    o.writer->open("TransDiffs_scatter.R");
    o.writer->write(RWriter::scatter(stats, ChrT, "", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
    o.writer->close();

    /*
     * Generating ROC plot
     */
    
    o.writer->open("TransDiffs_ROC.R");
    o.writer->write(RWriter::roc(stats.data.at(ChrT).seqs, stats.data.at(ChrT).qs));
    o.writer->close();

    /*
     * Generating LODR plot
     */

    o.writer->open("TransDiffs_ROC.R");
    o.writer->write(RWriter::roc(stats.data.at(ChrT).seqs, stats.data.at(ChrT).qs));
    o.writer->close();

    /*
     * Generating MA plot
     */

    o.writer->open("TransDiffs_ROC.R");
    o.writer->write(RWriter::roc(stats.data.at(ChrT).seqs, stats.data.at(ChrT).qs));
    o.writer->close();
    
}
