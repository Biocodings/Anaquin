/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "RnaQuin/r_diff.hpp"
#include "parsers/parser_adiff.hpp"

using namespace Anaquin;

extern Scripts PlotTROC();

typedef RDiff::Metrics Metrics;

std::vector<std::string> RDiff::classify(const std::vector<double> &qs, const std::vector<double> &folds, double qCut, double foldCut)
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

template <typename T> void classifyChrT(RDiff::Stats &stats, const T &t, const RDiff::Options &o)
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
                    known = match->concent(Mix_2) / match->concent(Mix_1);
                    
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
                    known = match->concent(Mix_2) / match->concent(Mix_1);
                    
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
    
    stats.elfs.push_back(log2(known));
}

template <typename T> void update(RDiff::Stats &stats, const T &x, const RDiff::Options &o)
{
    typedef DiffTest::Status Status;
    
    if (Standard::isSynthetic(x.cID))
    {
        stats.n_syn++;
        classifyChrT(stats, x, o);
    }
    else
    {
        stats.n_gen++;
        stats.elfs.push_back(NAN);
    }

    /*
     * The expected fold-change is done in the classifer because it can only happen with synthetic
     */

    stats.ids.push_back(x.id);

    stats.ps.push_back(x.p);
    stats.qs.push_back(x.q);
    stats.mlfs.push_back(x.logF);
    stats.ses.push_back(x.logFSE);
    stats.means.push_back(x.mean);
}

template <typename Functor> RDiff::Stats calculate(const RDiff::Options &o, Functor f)
{
    const auto &r = Standard::instance().r_trans;

    RDiff::Stats stats;

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
        case Metrics::Isoform: { stats.limit = r.absolute(stats.hist);     break; }
        default: { break; }
    }

    return stats;
}

RDiff::Stats RDiff::analyze(const std::vector<DiffTest> &tests, const Options &o)
{
    return calculate(o, [&](RDiff::Stats &stats)
    {
        for (auto &test : tests)
        {
            //update(stats, test, o);
        }
    });
}

RDiff::Stats RDiff::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](RDiff::Stats &stats)
    {
        ParserADiff::parse(file, [&](const ParserADiff::Data &x, const ParserProgress &)
        {
            update(stats, x, o);
        });
    });
}

void RDiff::report(const FileName &file, const Options &o)
{
    const auto m = std::map<Metrics, std::string>
    {
        { Metrics::Gene,    "genes"    },
        { Metrics::Isoform, "isoforms" },
    };

    const auto stats = RDiff::analyze(file, o);
    const auto units = m.at(o.metrs);
    
    o.info("Generating statistics");
    
    // Eg: DESeq2
    const auto shouldLODR = true;
    
    /*
     * Generating summary statistics
     */

    RDiff::generateSummary("RnaDiff_summary.stats", stats, o, units);

    /*
     * Generating differential results
     */

    RDiff::generateCSV("RnaDiff_quins.csv", stats, o);
    
    /*
     * Generating log-fold plot
     */
    
    o.generate("RnaDiff_fold.R");
    o.writer->open("RnaDiff_fold.R");
    o.writer->write(RWriter::createScatterNoLog("RnaDiff_quins.csv",
                                                "Fold Change",
                                                "Expected fold change (log2)",
                                                "Measured fold change (log2)",
                                                "Expected",
                                                "Measured", false));
    o.writer->close();

    /*
     * Generating ROC plot
     */
    
    o.generate("RnaDiff_ROC.R");
    o.writer->open("RnaDiff_ROC.R");
    o.writer->write(RWriter::createScript("RnaDiff_quins.csv", PlotTROC()));
    o.writer->close();

    /*
     * Generating LODR plot
     */
    
    if (shouldLODR)
    {
        RDiff::generateLODR("RnaDiff_LODR.R", "RnaDiff_quins.csv", o);
    }

    /*
     * Generating MA plot
     */
    
    //if (!o.counts.empty())
    //{
    //    RDiff::generateMA("RnaDiff_MA.R", "RnaDiff_counts.stats", o);
    //}
}