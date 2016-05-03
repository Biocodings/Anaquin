/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "data/pachter.hpp"
#include "TransQuin/t_kexpress.hpp"

using namespace Anaquin;

TKExpress::Stats TKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    TExpress::Options o_;
    
    o_.metrs = TExpress::Metrics::Isoform;
    o_.soft  = TExpress::Software::Kallisto;

    // Run quantification in Kallisto
    return TExpress::analyze(Pachter::externalQuant(o.index, file1, file2), o_);
}

void TKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto units = "isoforms";
    
    const auto files = std::vector<FileName> { file1, file2 };
    const auto stats = std::vector<TKExpress::Stats> { TKExpress::analyze(file1, file2, o) };
    
    /*
     * 1. Generating summary statistics (single or multiple samples)
     */
    
    TExpress::generateSummary("TransKExpress_summary.stats", files, stats, o, units);
    
    /*
     * 2. Generating detailed statistics
     */
    
    TExpress::generateCSV("TransKExpress_quins.csv", stats, o);
    
    /*
     * 3. Generating abundance vs abundance (single or multiple samples)
     */
    
    TExpress::generateRAbund("TransKExpress_express.R", "TransKExpress_quins.csv", stats, o);
    
    /*
     * 4. Generating major plot (but only if we have the isoforms...)
     */
    
    if (stats.size() >= 2)
    {
        TExpress::generateRSplice("TransKExpress_splice.R", "TransKExpress_quins.csv", o);
    }
}