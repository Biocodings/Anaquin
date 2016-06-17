/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "FusQuin/f_diff.hpp"
#include "FusQuin/FUSQuin.hpp"
#include "FusQuin/f_discover.hpp"
#include "parsers/parser_stab.hpp"
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotFFold();

FDiff::Stats FDiff::analyze(const FileName &normal, const FileName &fusion, const Options &o)
{
    FDiff::Stats stats;
    
    stats.normal = FNormal::analyze(normal);
    stats.fusion = FFusion::analyze(fusion);
    
    const auto &r = Standard::instance().r_fus;

    for (const auto &i : stats.normal)
    {
        for (const auto &j : stats.fusion)
        {
            if (FUSQuin::normToFusion(i.first) == j.first)
            {
                // Exptected relative expression
                const auto expected = r.findNormFus(i.first)->fold();
                
                // Measured relative expression
                const auto measured = static_cast<double>(i.second.y) / j.second.y;
  
                stats.add(i.first + " - " + j.first, expected, measured);
            }
        }
    }
    
    stats.hist = stats.normal.hist;
    stats.hist.insert(stats.fusion.hist.begin(), stats.fusion.hist.end());

    if (stats.normal.limit.abund < stats.fusion.limit.abund)
    {
        stats.limit = stats.normal.limit;
    }
    else
    {
        stats.limit = stats.fusion.limit;
    }
    
    return stats;
}

void FDiff::report(const FileName &normal, const FileName &fusion, const Options &o)
{
    const auto stats = FDiff::analyze(normal, fusion, o);

    o.info("Generating statistics");
    
    /*
     * Generating FusDiff_summary.stats
     */
    
    o.writer->open("FusDiff_summary.stats");
    o.writer->write(StatsWriter::linearSummary(normal + " & " + fusion, o.rAnnot, stats, stats, stats.hist, "sequins"));
    o.writer->close();
    
    /*
     * Generating FusDiff_quins.stats
     */
    
    o.info("Generating FusDiff_quins.csv");
    o.writer->open("FusDiff_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating FusDiff_fold.R
     */
    
    o.info("Generating FusDiff_fold.R");
    o.writer->open("FusDiff_fold.R");
    o.writer->write(RWriter::createScript("FusDiff_quins.csv", PlotFFold()));
    o.writer->close();
}