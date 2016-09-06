#include "tools/script.hpp"
#include "RnaQuin/r_kexpress.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts KExpressScript();

RKExpress::Stats RKExpress::analyze(const FileName &index, const FileName &p1, const FileName &p2, const Options &o)
{
    // Where the analysis should be saved
    const auto output = tmpFile();
    
    // The script to run
    const auto script = KExpressScript();
    
    // Eg: python kexpress.py RNAQuin ARN004.v032.index /tmp/kallisto LRN087.1_val_1.fq LRN087.2_val_2.fq
    Script::run(script, "python", "RnaQuin " + index + " " + output + " " + p1 + " " + p2);
    
    RExpress::Options ro;
    
    ro.writer = o.writer;
    ro.metrs  = RExpress::Metrics::Gene;
    ro.format = RExpress::Format::Kallisto;
    
    const auto abund = output + "/abundance.tsv";
    
    RKExpress::Stats stats(RExpress::analyze(abund, ro));
    stats.abund = abund;

    return stats;
}

void RKExpress::report(const FileName &index, const FileName &p1, const FileName &p2, const Options &o)
{
    const auto stats = RKExpress::analyze(index, p1, p2, o);
    
    RExpress::Options ro;
    ro.writer = o.writer;

    /*
     * Generating RnaKExpress_summary.stats
     */
    
    RExpress::generateSummary("RnaKExpress_summary.stats", std::vector<FileName> { stats.abund }, std::vector<RExpress::Stats> { stats.stats }, ro, "genes");
}