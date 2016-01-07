/*
 * Copyright (C) 2015 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "trans/t_diffs.hpp"
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

typedef TDiffs::Metrics  Metrics;
typedef TDiffs::Software Software;

template <typename T> void update(TDiffs::Stats &stats, const T &t, const GenericID &id, bool isoform, const TDiffs::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    
    if (t.cID != ChrT)
    {
        stats.n_expT++;
    }
    else
    {
        stats.n_chrT++;
    }
    
    // The known and observed fold-change
    Fold known = NAN;
    
    // It's NAN if the sequin defined in reference but not in mixture
    Fold measured = NAN;
    
    /*
     * Differential expression at the gene level
     */
    
    auto g = [&](const GeneID &id, double fpkm_1, double fpkm_2)
    {
        // The known and observed fold-change
        Fold known = NAN;
        
        // It's NAN if the sequin defined in reference but not in mixture
        Fold measured = NAN;
        
        const auto *g = r.findGene(t.cID, id);
        
        if (g)
        {
            // Calculate the known fold-change between B and A
            known = (g->abund(Mix_2) / g->abund(Mix_1));
        }
        
        if (g && !isnan(fpkm_1) && !isnan(fpkm_2) && fpkm_1 && fpkm_2)
        {
            stats.h.at(id)++;
            
            // Measured fold-change between the two mixtures
            measured = fpkm_2 / fpkm_1;
        }
        
        stats.data[t.cID].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
    };

    switch (o.metrs)
    {
        case Metrics::Gene:
        {
            if ((t.status != NoTest) && stats.h.count(t.id))
            {
                g(t.id, t.fpkm_1, t.fpkm_2);
            }
            
            break;
        }
            
        case Metrics::Isoform:
        {
            if ((t.status == NoTest) || !stats.h.count(id))
            {
                return;
            }
            
            const auto *seq = r.match(id);
            
            if (seq)
            {
                // Known fold-change between the two mixtures
                known = seq->abund(Mix_2) / seq->abund(Mix_1);
            }
            
            if ((t.status != NoTest) && t.fpkm_1 && t.fpkm_2)
            {
                stats.h.at(id)++;
                
                // Measured fold-change between the two mixtures
                measured = t.fpkm_2 / t.fpkm_1;
            }
            
            stats.data[t.cID].add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            
            break;
        }
    }
}

template <typename Functor> TDiffs::Stats calculate(const TDiffs::Options &o, Functor f)
{
    TDiffs::Stats stats;

    const auto &r = Standard::instance().r_trans;
    const auto cIDs = r.chromoIDs();
    
    std::for_each(cIDs.begin(), cIDs.end(), [&](const ChromoID &cID)
    {
        stats.data[cID];
    });

    const auto isoform = (o.metrs == Metrics::Isoform);
    o.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    // Construct for a histogram at the appropriate level
    stats.h = isoform ? r.hist() : r.geneHist(ChrT);
    
    o.info("Parsing tracking file");

    f(stats);
    
    /*
     * Calculating detection limit
     */
    
    o.info("Calculating detection limit");
    
    stats.s = isoform ? r.limit(stats.h) : r.limitGene(stats.h);
    
    return stats;
}

TDiffs::Stats TDiffs::analyze(const std::vector<DiffTest> &tests, const Options &o)
{
    return calculate(o, [&](TDiffs::Stats &stats)
    {
        for (auto &test : tests)
        {
            update(stats, test, test.id, o.metrs == Metrics::Isoform, o);
        }
    });
}

TDiffs::Stats TDiffs::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](TDiffs::Stats &stats)
    {
        switch (o.soft)
        {
            case Software::Cuffdiffs:
            {
                const auto isIsoform = (o.metrs == Metrics::Isoform);
                
                ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
                {
                    update(stats, t, isIsoform ? t.testID : t.id, isIsoform, o);
                });

                break;
            }
        }
    });
}

void TDiffs::report(const FileName &file, const Options &o)
{
    const auto stats = TDiffs::analyze(file, o);
    const auto units = (o.metrs == Metrics::Isoform) ? "isoforms" : "genes";
    
    o.info("Generating statistics");
    
    for (const auto &i : stats.data)
    {
        /*
         * Generating summary statistics
         */
        
        o.writer->open("TransDiffs_summary.stats");
        o.writer->write(RWriter::linear(file, stats, i.first, units));
        o.writer->close();
        
        /*
         * Generating scatter plot
         */
        
        o.writer->open("TransDiffs_scatter.R");
        o.writer->write(RWriter::scatter(stats, i.first, "", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
        o.writer->close();
    }
}

void TDiffs::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = analyze(files, o);
    //const auto units = (o.level == Isoform) ? "isoforms" : "genes";

    /*
     * Generating summary statistics for each replicate
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        //const auto file = (boost::format("TransDiff_%1%_summary.stats") % files[i]).str();
        ////AnalyzeReporter::linear(file, files[i], stats[i], units, o.writer);
    }
    
    /*
     * Generating scatter plots for each replicate
     */

    for (auto i = 0; i < files.size(); i++)
    {
        //const auto file = (boost::format("TransDiff_%1%_summary.stats") % files[i]).str();
        ////AnalyzeReporter::scatter(stats[i], "", file, "Expected fold change of mixture A and B", "Measured fold change of mixture A and B", "Expected log2 fold change of mixture A and B", "Expected log2 fold change of mixture A and B", o.writer);
    }
    
    /*
     * Generating summary statistics for all replicates
     */
    
}