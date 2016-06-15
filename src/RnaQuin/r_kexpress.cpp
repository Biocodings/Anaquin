/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "data/pachter.hpp"
#include "RnaQuin/r_kexpress.hpp"

using namespace Anaquin;

TKExpress::Stats TKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    TExpress::Options o_;
    
    o_.metrs = TExpress::Metrics::Isoform;
    //o_.soft  = TExpress::Software::Kallisto;

    // Run quantification in Kallisto
    return TExpress::analyze(Pachter::externalQuant(o.index, file1, file2), o_);
}

void TKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    //const auto units = "isoforms";
    
    const auto files = std::vector<FileName> { file1, file2 };
    const auto stats = std::vector<TKExpress::Stats> { TKExpress::analyze(file1, file2, o) };
    
    /*
     * Generating RnaKExpress_summary.stats
     */
    
    //TExpress::generateSummary("RnaKExpress_summary.stats", files, stats, o, units);
    
    /*
     * Generating RnaKExpress_quins.csv
     */
    
    TExpress::generateCSV("RnaKExpress_quins.csv", stats, o);
    
    /*
     * Generating RnaKExpress_express.R
     */
    
    //TExpress::generateR("RnaKExpress_express.R", "RnaKExpress_quins.csv", stats, o);
    
    /*
     * Generating RnaKExpress_splice.R
     */
    
    if (stats.size() >= 2)
    {
        TExpress::generateRSplice("RnaKExpress_splice.R", "RnaKExpress_quins.csv", o);
    }
}