/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "trans/t_diffs.hpp"
#include "parsers/parser_edgeR.hpp"
#include "parsers/parser_DESeq2.hpp"
#include "parsers/parser_cdiffs.hpp"

using namespace Anaquin;

typedef TDiffs::Level    Level;
typedef TDiffs::Software Software;
typedef DiffTest::Status Status;

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
    
    const auto &id = t.id;
    
    const auto &r = Standard::instance().r_trans;
    
    // Known fold change
    Fold known = NAN;
    
    // It's NAN if the sequin defined in reference but not in mixture
    Fold measured = NAN;

    switch (o.lvl)
    {
        case Level::Gene:
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
            
        case Level::Isoform:
        {
            if (t.status == Status::Tested && !stats.hist.count(id))
            {
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
                    measured = t.logF;

                    // Turn it back to the original scale
                    measured = std::pow(2, measured);
                }
                
                stats.data[ChrT].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            }

            break;
        }
            
        case Level::Exon:
        {
            throw "Not Implemented";
        }
    }
}

template <typename T> void classifyEndo(TDiffs::Stats &, const T &, const TDiffs::Options &)
{
    /*
     * Obviously, we can't compare with the expected here...
     */
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

    stats.data[t.cID].ps.push_back(t.p);
    stats.data[t.cID].qs.push_back(t.q);
    stats.data[t.cID].ids.push_back(t.id);
    stats.data[t.cID].logFs.push_back(t.logF);
}

template <typename Functor> TDiffs::Stats calculate(const TDiffs::Options &o, Functor f)
{
    TDiffs::Stats stats;

    const auto &r = Standard::instance().r_trans;
    const auto cIDs = r.chromoIDs();

    stats.data[Endo];
    stats.data[ChrT];

    switch (o.lvl)
    {
        case Level::Gene:    { stats.hist = r.geneHist(ChrT); break; }
        case Level::Isoform: { stats.hist = r.hist();         break; }
        case Level::Exon:    { throw "Not Implemented"; }
    }

    assert(!stats.hist.empty());
    
    o.info("Parsing input files");
    f(stats);
    
    /*
     * Calculating detection limit
     */
    
    o.info("Calculating detection limit");
    
    switch (o.lvl)
    {
        case Level::Gene:    { stats.limit = r.limitGene(stats.hist); break; }
        case Level::Isoform: { stats.limit = r.limit(r.hist());       break; }
        case Level::Exon:
        {
            throw "Not Implemented";
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
        switch (o.soft)
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
                ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
                {
                    update(stats, t, o);
                });

                break;
            }
        }
    });
}

void TDiffs::report(const FileName &file, const Options &o)
{
    const auto stats = TDiffs::analyze(file, o);
    
    const auto m = std::map<TDiffs::Level, std::string>
    {
        { TDiffs::Level::Gene, "gene"    },
        { TDiffs::Level::Gene, "exon"    },
        { TDiffs::Level::Gene, "isoform" },
    };
    
    const auto units = m.at(o.lvl);
    
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
     * Generating scatter plot for the log-fold changes
     */
    
    o.writer->open("TransDiffs_scatter.R");
    o.writer->write(RWriter::scatter(stats, ChrT, "", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
    o.writer->close();

    /*
     * Generating ROC plot
     */
    
    o.writer->open("TransDiffs_ROC.R");
    o.writer->write(RWriter::roc(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps));
    o.writer->close();

    /*
     * Generating LODR plot
     */
    
 //   o.writer->open("TransDiffs_LODR.R");
   // o.writer->write(RWriter::lodr(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps));
   // o.writer->close();

    
    /*
     * Generating MA plot
     */

//    o.writer->open("TransDiffs_MA.R");
  //  o.writer->write(RWriter::roc(stats.data.at(ChrT).ids, stats.data.at(ChrT).qs));
    //o.writer->close();
}
