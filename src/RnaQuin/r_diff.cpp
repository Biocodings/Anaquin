/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "RnaQuin/r_diff.hpp"
#include "parsers/parser_diff.hpp"

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

template <typename T> void classifySyn(RDiff::Stats &stats, const T &t, const RDiff::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    assert(Standard::isSynthetic(t.cID));
    
    SequinID id;
    
    Fold exp = NAN;
    Fold obs = NAN;

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
                
                id = t.id;
                
                // Calculate the known fold-change between B and A
                exp = exp_2 / exp_1;
                
                // Measured fold-change between the two mixtures
                obs = t.logF;
                
                // Turn it back to the original scale
                obs = std::pow(2, obs);
            }

            break;
        }

        case Metrics::Isoform:
        {
            const auto match = r.match(t.id);
            
            if (match)
            {
                stats.hist.at(t.cID).at(t.id)++;
                
                id = t.id;
                
                // Known fold-change between the two mixtures
                exp = match->concent(Mix_2) / match->concent(Mix_1);
                
                // Measured fold-change between the two mixtures
                obs = t.logF;
                
                // Turn it back to the original scale
                obs = std::pow(2, obs);
            }
            
            break;
        }
    }
    
    if (!id.empty())
    {
        // This is not on the log scale, so it can't be non-positive
        assert(exp > 0);

        stats.add(id, !isnan(exp) ? exp : NAN, !isnan(obs) ? obs : NAN);
        
        if (isnan(stats.limit.abund) || exp < stats.limit.abund)
        {
            stats.limit.id = id;
            stats.limit.abund = exp;
        }
        
        stats.elfs.push_back(log2(exp));
    }
    else
    {
        o.warn(t.id + " not found");
    }
}

template <typename T> void update(RDiff::Stats &stats, const T &x, const RDiff::Options &o)
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
        stats.elfs.push_back(NAN);
    }

    /*
     * The expected fold-change is done in the classifer because it can only happen with synthetic
     */

    stats.ps.push_back(x.p);
    stats.qs.push_back(x.q);
    stats.ids.push_back(x.id);
    stats.cIDs.push_back(x.cID);
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
        case Metrics::Gene:    { stats.hist = r.histGene(); break; }
        case Metrics::Isoform: { stats.hist = r.histIsof(); break; }
    }

    assert(!stats.hist.empty());

    o.info("Parsing input files");
    f(stats);
    
    return stats;
}

RDiff::Stats RDiff::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](RDiff::Stats &stats)
    {
        ParserDiff::parse(file, [&](const ParserDiff::Data &x, const ParserProgress &)
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
    const auto shouldLOD = true;
    
    /*
     * Generating RnaFoldChange_summary.stats
     */

    RDiff::generateSummary("RnaFoldChange_summary.stats", stats, o, units);

    /*
     * Generating RnaFoldChange_sequins.csv
     */

    RDiff::generateCSV("RnaFoldChange_sequins.csv", stats, o);
    
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
        RDiff::generateLODR("RnaFoldChange_LODR.R", "RnaFoldChange_sequins.csv", o);
    }

    /*
     * Generating MA plot
     */
    
    //if (!o.counts.empty())
    //{
    //    RDiff::generateMA("RnaFoldChange_MA.R", "RnaFoldChange_counts.stats", o);
    //}
}