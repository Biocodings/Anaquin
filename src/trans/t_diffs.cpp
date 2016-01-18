/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "trans/t_diffs.hpp"
#include "data/experiment.hpp"
#include "parsers/parser_edgeR.hpp"
#include "parsers/parser_DESeq2.hpp"
#include "parsers/parser_cdiffs.hpp"
#include "parsers/parser_HTSeqCount.hpp"

using namespace Anaquin;

typedef TDiffs::Level    Level;
typedef DiffTest::Status Status;
typedef TDiffs::DiffSoft Software;

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

static void counts(TDiffs::Stats &stats, const TDiffs::Options &o)
{
    // Create an empty count table
    stats.counts = std::shared_ptr<CountTable>(new CountTable(o.exp->countTable()));

    const auto &names = stats.counts->names();

    assert(o.counts.size() == names.size());
    
    for (auto i = 0; i < o.counts.size(); i++)
    {
        switch (o.cSoft)
        {
            case TDiffs::CountSoft::HTSeqCount:
            {
                const auto &name = names.at(i);
                
                ParserHTSeqCount::parse(Reader(o.counts[i]), [&](const ParserHTSeqCount::Sample &s, const ParserProgress &)
                {
                    if (!i)
                    {
                        stats.counts->addFeature(s.id);
                    }
                    
                    stats.counts->addCount(name, s.count);
                });

                break;
            }
        }
    }
    
    
    
    
//    
//    
//    /*
//     * Create a reader for each replicate, we'll read them simultaneously
//     */
//    
//    std::vector<Reader> rs;
//    
//    for (auto i = 0; i < o.counts.size(); i++)
  //  {
    //    rs.push_back(Reader(o.counts.at(i)));
    //}
//    
//    // One for each condition, sorted by the factor level
//    stats.avgs.resize(o.exp->countConds());
//    
//    /*
//     * Before we begin, we should cache the factors for each condition
//     */
//    
//    auto conds = std::map<Experiment::Factor, std::vector<std::size_t>>
//    {
//        { 0, o.exp->cond(0) }, { 1, o.exp->cond(1) }
//    };
//    
//    // This is differential analysis, thus we must always have two conditions...
//    assert(conds.size() == 2);
//    
//    /*
//     * Read the count files and calculate their arithemetic averages
//     */
//    
//    switch (o.cSoft)
//    {
//        case TDiffs::CountSoft::HTSeqCount:
//        {
//            
//            
//            
//            
////            ParserHTSeqCount::parse(rs, [&](const ParserHTSeqCount::Samples &s, const ParserProgress &)
////            {
////                /*
////                 * We shouldn't assume the orders of the conditions. For example, we could be given:
////                 *
////                 *      A1,A2,A3,B1,B2,B2 or B1,B3,B2,A3,A2,A1
////                 *
////                 * for the same experiment.
////                 */
////                
////                for (auto i = 0; i < conds.size(); i++)
////                {
////                    std::vector<unsigned> counts;
////                    
////                    for (auto j = 0; j < conds[i].size(); j++)
////                    {
////                        counts.push_back(s.counts[conds[i][j]]);
////                    }
////                    
////                    stats.avgs[i][s.id] = SS::mean(counts);
////                }
////            });
//            
//            break;
//        }
 //   }
}

template <typename Functor> TDiffs::Stats calculate(const TDiffs::Options &o, Functor f)
{
    TDiffs::Stats stats;

    const auto &r = Standard::instance().r_trans;

    stats.data[Endo];
    stats.data[ChrT];

    switch (o.lvl)
    {
        case Level::Gene:    { stats.hist = r.geneHist(ChrT); break; }
        case Level::Isoform: { stats.hist = r.hist();         break; }
        case Level::Exon:    { throw "Not Implemented"; }
    }

    assert(!stats.hist.empty());

    /*
     * Calculate average counts for each condition. This is optional but needed for MA plot and LODR plot.
     */
    
    if (!o.counts.empty())
    {
        counts(stats, o);
    }
    
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
        switch (o.dSoft)
        {
            case DiffSoft::DESeq2:
            {
                ParserDESeq2::parse(file, [&](const DiffTest &t, const ParserProgress &)
                {
                    update(stats, t, o);
                });

                break;                
            }
                
            case DiffSoft::edgeR:
            {
                break;
            }

            case DiffSoft::Cuffdiff:
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

static void writeCounts(const FileName &file, const TDiffs::Stats &stats, const TDiffs::Options &o)
{
    o.writer->open(file);

    const auto &names = stats.counts->names();
    
    for (const auto &name : names)
    {
        o.writer->write("," + name, false);
    }

    o.writer->write("\n", false);
    
    const auto &ids = stats.counts->ids();

    for (auto i = 0; i < ids.size(); i++)
    {
        o.writer->write(ids[i], false);

        for (const auto &name : names)
        {
            o.writer->write("," + std::to_string(stats.counts->counts(name).at(i)), false);
        }

        o.writer->write("\n", false);
    }
    
    o.writer->close();
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
    o.writer->write(RWriter::createROC(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps));
    o.writer->close();

    /*
     * Generating count table
     */

    writeCounts("TransDiffs_counts.csv", stats, o);

    /*
     * Generating MA plot
     */

    o.writer->open("TransDiffs_MA.R");
    o.writer->write(RWriter::createMA("TransDiffs_counts.csv", units));
    o.writer->close();
    
    /*
     * Generating LODR plot
     */
    
    //   o.writer->open("TransDiffs_LODR.R");
    // o.writer->write(RWriter::lodr(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps));
    // o.writer->close();
    
    
}
