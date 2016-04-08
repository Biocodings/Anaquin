/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "stats/kallisto.hpp"
#include "TransQuin/t_express.hpp"
#include "TransQuin/t_kexpress.hpp"

using namespace Anaquin;

TKExpress::Stats TKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_trans;
    
    TKExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();
    
    stats.n_endo = NAN;
    
    // Run quantification in Kallisto
    const auto abundFile = Kallisto::quant(o.index, file1, file2);
    
    
    

    
    
    return stats;
}

void TKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto stats = TKExpress::analyze(file1, file2, o);
    
//    o.info("Generating statistics");
//    
//    const auto units = "isoforms";
//    
//    /*
//     * 1. Generating summary statistics (single or multiple samples)
//     */
//    
//    o.info("Generating TransKExpress_summary.stats");
//    o.writer->open("TransKExpress_summary.stats");
//    
//    if (files.size() == 1)
//    {
//        o.writer->write(TExpress::singleSummary(stats[0], files[0], units, o));
//    }
//    else
//    {
//        o.writer->write(TExpress::multipleSummary(files, stats, units, o));
//    }
//
//    o.writer->close();
//    
//    /*
//     * 2. Generating detailed statistics for the sequins
//     */
//    
//    o.info("Generating TransKExpress_quins.csv");
//    o.writer->open("TransKExpress_quins.csv");
//    
//    if (files.size() == 1)
//    {
//        o.writer->write(StatsWriter::writeCSV(stats[0], "EAbund", "MAbund"));
//    }
//    else
//    {
//        o.writer->write(TExpress::multipleCSV(stats));
//    }
//    
//    o.writer->close();
//    
//    /*
//     * 3. Generating abundance vs abundance (single or multiple samples)
//     */
//    
//    o.info("Generating TransKExpress_abundAbund.R");
//    o.writer->open("TransKExpress_abundAbund.R");
//    
//    if (files.size() == 1)
//    {
//        o.writer->write(RWriter::createScript("TransKExpress_quins.csv", PlotTAbundAbund()));
//    }
//    else
//    {
//        o.writer->write(RWriter::createScript("TransKExpress_quins.csv", PlotRAbundAbund()));
//    }
//    
//    o.writer->close();
//    
//    /*
//     * 4. Generating major plot (but only if we have the isoforms...)
//     */
//
//    if (files.size() >= 2)
//    {
//        o.writer->open("TransKExpress_major.R");
//        o.writer->write(RWriter::createScript("TransKExpress_quins.csv", PlotMajor()));
//        o.writer->close();
//    }
}