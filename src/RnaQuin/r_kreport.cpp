#include "tools/script.hpp"
#include "RnaQuin/r_kreport.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

// Shared with other modules
extern Path __output__;

// Defined in resources.cpp
extern Scripts KExpressScript();

RKReport::Stats RKReport::analyze(const FileName &index, const FileName &p1, const FileName &p2, const Options &o)
{
    // Where the analysis should be saved
    const std::string output = tmpFile();
    
    // The script to run
    const auto script = KExpressScript();
    
    // Eg: python kexpress.py RNAQuin ARN004.v032.index /tmp/kallisto LRN087.1_val_1.fq LRN087.2_val_2.fq
    Script::run(script, "python", "RnaQuin " + index + " " + output + " " + p1 + " " + p2, ".py"    );
    
    RExpress::Options ro;
    
    ro.writer = o.writer;
    ro.metrs  = RExpress::Metrics::Gene;
    ro.format = RExpress::Format::Kallisto;
    
    const auto abund = output + "/abundance.tsv";
    
    RKReport::Stats stats(RExpress::analyze(abund, ro));
    stats.abund  = abund;
    stats.output = output;

    return stats;
}

void RKReport::report(const FileName &index, const FileName &p1, const FileName &p2, const Options &o)
{
    const auto stats = RKReport::analyze(index, p1, p2, o);
    
    o.info("Temporary: " + stats.output);
    
    RExpress::Options ro;

    // Make sure to reuse the temporary directory
    ro.writer = std::shared_ptr<FileWriter>(new FileWriter(stats.output));
    
    // We're doing gene expression
    ro.metrs = RExpress::Metrics::Gene;
    
    ro.logger = o.logger;
    
    const auto vStats = std::vector<RExpress::Stats> { stats.stats };
    
    __output__ = stats.output;
    
    /*
     * Generating outputs to the temporatory directory, from which we can create our report.
     */
    
    RExpress::generateCSV("RnaKReport_sequins.csv", vStats, ro);
    RExpress::generateSummary("RnaKReport_summary.stats", std::vector<FileName> { stats.abund }, vStats, ro, "genes");
    RExpress::generateR("RnaKReport_linear.R", "RnaKReport_sequins.csv", vStats, ro);

    /*
     * Create a PDF report based on those files
     */

    Script::report("RnaQuin", "RnaKReport_report.pdf", stats.output, o);
}