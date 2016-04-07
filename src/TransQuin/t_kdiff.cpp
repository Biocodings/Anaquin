/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "TransQuin/t_kdiff.hpp"
#include "parsers/parser_sleuth.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotFold();

// Defined in resources.cpp
extern Scripts PlotLODR();

// Defined in resources.cpp
extern Scripts PlotTROC();

typedef DiffTest::Status Status;
typedef TKDiff::Software Software;

template <typename T> void match(TKDiff::Stats &stats, const T &t, const TKDiff::Options &o)
{
    const auto &id = t.id;
    const auto &r  = Standard::instance().r_trans;
    
    // Known fold change
    Fold known = NAN;
    
    // It's NAN if the sequin defined in reference but not in mixture
    Fold measured = NAN;
    
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
}

TKDiff::Stats TKDiff::analyze(const FileName &file, const Options &o)
{
    const auto &r  = Standard::instance().r_trans;
    
    TKDiff::Stats stats;
    stats.hist = r.hist();
    
    ParserSleuth::parse(file, [&](const ParserSleuth::Data &data, const ParserProgress &)
    {
        match(stats, data, o);
    });
    
    return stats;
}

void TKDiff::report(const FileName &file, const Options &o)
{
    const auto stats = TKDiff::analyze(file, o);
    //const auto units = "isoform";

    o.info("Generating statistics");
    
    /*
     * 1. Generating summary statistics
     */
    
    o.info("Generating TransKDiff_summary.stats");
    o.writer->open("TransKDiff_summary.stats");
    o.writer->write(StatsWriter::linearSummary(file, o.rChrT, stats, stats.hist));
    o.writer->close();

    /*
     * 2. Generating differential results
     */
    
    o.info("Generating TransKDiff_quins.csv");
    o.writer->open("TransKDiff_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats, "ELFold", "MLFold", true));
    o.writer->close();
    
    /*
     * 3. Generating log-fold plot
     */
    
    o.info("Generating TransKDiff_fold.R");
    o.writer->open("TransKDiff_fold.R");
    o.writer->write(RWriter::createScript("TransKDiff_quins.csv", PlotFold()));
    o.writer->close();
}
