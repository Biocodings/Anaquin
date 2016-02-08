#include "fusion/FUSQuin.hpp"
#include "fusion/f_express.hpp"

using namespace Anaquin;

FExpress::Stats FExpress::analyze(const FileName &file, const FDiscover::Options &o)
{
    FExpress::Stats stats;
    
    stats.data[ChrT];
    
    return stats;
}

void FExpress::report(const FileName &file, const FDiscover::Options &o)
{
    const auto stats = analyze(file, o);

//    o.logInfo("Calculating limit of sensitivity");
//    stats.ss = Standard::instance().r_fus.limit(stats.chrT->h);
//
//    /*
//     * Generating summary statistics
//     */
//
//    o.info("Generating summary statistics");
//    //AnalyzeReporter::linear("FusionExpress_summary.stats", file, stats, "fusions", o.writer);
//
//    /*
//     * Generating Bioconductor
//     */
//
//    o.info("Generating Bioconductor");
//    //AnalyzeReporter::scatter(stats, "", "FusionExpress", "Expected concentration (attomol/ul)", "Measured coverage (reads)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 reads)", o.writer);
//    
}