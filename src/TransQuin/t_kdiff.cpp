/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "TransQuin/t_kdiff.hpp"
#include "parsers/parser_sleuth.hpp"

using namespace Anaquin;

typedef DiffTest::Status Status;
typedef TKDiff::Software Software;

TKDiff::Stats TKDiff::analyze(const FileName &file, const Options &o)
{
    TKDiff::Stats stats;
    
    ParserSleuth::parse(file, [&](const ParserSleuth::Data &d, const ParserProgress &)
    {
                            
    });
    
    return stats;
}

void TKDiff::report(const FileName &file, const Options &o)
{
//    const auto stats = TKDiff::analyze(file, o);
//
//    const auto units = m.at(o.metrs);
//    
//    o.info("Generating statistics");
//    
//    /*
//     * There's no need to write summary for each replicate because differential analysis has already incorporated them.
//     */
//    
//    /*
//     *  - Summary statistics
//     *  - Sequin statistics
//     *  - Scatter plot
//     *  - ROC plot
//     *  - MA plot
//     *  - LODR Plot
//     */
//
//    /*
//     * Generating summary statistics
//     */
//    
//    o.writer->open("TransDiff_summary.stats");
//    //o.writer->write(StatsWriter::linear(file, stats, ChrT, units));
//    o.writer->close();
//    
//    /*
//     * Generating scatter plot for the log-fold changes
//     */
//    
//    o.writer->open("TransDiff_scatter.R");
//    o.writer->write(RWriter::scatter(stats, ChrT, "????", "TransDiff", "Expected fold change", "Measured fold change", "Expected log2 fold change", "Measured log2 fold change"));
//    o.writer->close();
//
//    /*
//     * Generating ROC plot
//     */
//    
//    o.writer->open("TransDiff_ROC.R");
//    //o.writer->write(RWriter::createROC_T(stats.data.at(ChrT).ids, stats.data.at(ChrT).ps, units));
//    o.writer->close();
//
//    /*
//     * Generating differential results (CSV)
//     */
//    
//    writeDifferent("TransDiff_diffs.csv", stats, o);
//
//    /*
//     * Generating MA plot
//     */
//    
//    o.writer->open("TransDiff_MA.R");
//    o.writer->write(RWriter::createMA("TransDiff_diffs.csv", units));
//    o.writer->close();
//    
//    /*
//     * Generating LODR plot
//     */
//    
//    o.writer->open("TransDiff_LODR.R");
//    o.writer->write(RWriter::createLODR_T("TransDiff_diffs.csv"));
//    o.writer->close();
}
