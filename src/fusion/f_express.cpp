#include "fusion/f_express.hpp"
#include "fusion/f_classify.hpp"

using namespace Anaquin;

FExpress::Stats FExpress::report(const FileName &file, const FDiscover::Options &o)
{
    auto stats = FClassify::analyze<FExpress::Options, FExpress::Stats>(file, true, o);

    o.logInfo("Calculating limit of sensitivity");
    stats.ss = Standard::instance().r_fus.limit(stats.h);

    /*
     * Generate summary statistics
     */

    o.info("Generating summary statistics");
    AnalyzeReporter::linear("FusionExpress_summary.stats", stats, "fusions", o.writer);

    /*
     * Generate an R script
     */

    o.info("Generating R script");
    AnalyzeReporter::scatter(stats, "", "FusionExpress", "Expected abudnance (log2 attomol/ul)", "Measured coverage (log2 reads)", o.writer);
    
    return stats;
}