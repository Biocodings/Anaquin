/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "TransQuin/t_kdiff.hpp"

using namespace Anaquin;

typedef DiffTest::Status Status;

TKDiff::Stats TKDiff::analyze(const std::vector<FileName> &files, const Options &o)
{
    const auto &r  = Standard::instance().r_trans;

    TKDiff::Stats stats;

    // Initialize the distribution for each sequin
    stats.hist = r.hist();
    
    return stats;
}

void TKDiff::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TKDiff::analyze(files, o);
    //const auto units = "isoform";

//    o.info("Generating statistics");
//    
//    /*
//     * 1. Generating summary statistics
//     */
//    
//    o.info("Generating TransKDiff_summary.stats");
//    o.writer->open("TransKDiff_summary.stats");
//    o.writer->write(StatsWriter::linearSummary(file, o.rChrT, stats, stats.hist));
//    o.writer->close();
//
//    /*
//     * 2. Generating differential results
//     */
//    
//    o.info("Generating TransKDiff_quins.csv");
//    o.writer->open("TransKDiff_quins.csv");
//    o.writer->write(StatsWriter::writeCSV(stats, "ELFold", "MLFold", true));
//    o.writer->close();
//    
//    /*
//     * 3. Generating log-fold plot
//     */
//    
//    o.info("Generating TransKDiff_fold.R");
//    o.writer->open("TransKDiff_fold.R");
//    o.writer->write(RWriter::createScript("TransKDiff_quins.csv", PlotFold()));
//    o.writer->close();
}
