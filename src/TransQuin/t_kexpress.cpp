/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "stats/kallisto.hpp"
#include "TransQuin/t_kexpress.hpp"

using namespace Anaquin;

TKExpress::Stats TKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_trans;
    
    TKExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();
    
    stats.n_endo = NAN;
    
    TExpress::Options o_;
    
    o_.soft = TExpress::Software::Kallisto;

    // Run quantification in Kallisto
    return TExpress::analyze(Kallisto::quant(o.index, file1, file2), o_);
}

void TKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto units = "isoforms";
    
    const auto files = std::vector<FileName> { file1, file2 };
    const auto stats = std::vector<TKExpress::Stats> { TKExpress::analyze(file1, file2, o) };
    
    /*
     * 1. Generating summary statistics (single or multiple samples)
     */
    
    TExpress::generateSummary("TransExpress_summary.stats", files, stats, o, units);
    
    /*
     * 2. Generating detailed statistics for the sequins
     */
    
    TExpress::generateCSV("TransExpress_quins.csv", stats, o);
    
    /*
     * 3. Generating abundance vs abundance (single or multiple samples)
     */
    
    TExpress::generateRAbund("TransExpress_abundAbund.R", "TransExpress_quins.csv", stats, o);
    
    /*
     * 4. Generating major plot (but only if we have the isoforms...)
     */
    
    TExpress::generateRMajor("TransExpress_major.R", "TransExpress_quins.csv", stats, o);
}