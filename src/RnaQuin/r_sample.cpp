#include "tools/sample.hpp"
#include "RnaQuin/r_sample.hpp"
#include "writers/writer_sam.hpp"

using namespace Anaquin;

RSample::Stats RSample::stats(const FileName &file, const Options &o)
{
    RSample::Stats stats;
    
    o.info(file);
    
    switch (o.meth)
    {
        case Method::Prop: { stats.p = o.p; break; }
    }

    assert(stats.p > 0 && stats.p < 1.0);
    o.info("Sampling proportion: " + std::to_string(stats.p));

    // Perform subsampling
    Sampler::subsample(file, stats.p, o);

    return stats;
}

static void generateSummary(const FileName &file, const RSample::Stats &stats, const RSample::Options &o)
{
    o.generate(file);
    
    auto meth2Str = [&]()
    {
        switch (o.meth)
        {
            case RSample::Method::Prop:  { return std::to_string(o.p); }
        }
    };
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Proportion of alignments mapped to the synthetic and genome\n"
                         "   ***\n\n"
                         "   Unmapped:  %2% reads\n"
                         "   Synthetic: %3% reads\n"
                         "   Genome:    %4% reads\n\n"
                         "   ***\n"
                         "   *** Reference annotation\n"
                         "   ***\n\n"
                         "   File: %5%\n\n"
                         "   Reference: %6% sequins\n\n"
                         "   ***                                               \n"
                         "   *** How the coverage is calculated?               \n"
                         "   ***                                               \n"
                         "   *** Please refer to the documentation for details \n"
                         "   ***                                               \n\n"
                         "   Method: %7%\n\n"
                         "   ***                                     \n"
                         "   ***    Statistics before subsampling    \n"
                         "   ***                                     \n\n"
                         "   Aligns: %8%\n"
                         "   Coverage (Synthetic): %9%\n"
                         "   Coverage (Genome):    %10%\n\n"
                         "   ***                                     \n"
                         "   ***    Statistics after subsampling     \n"
                         "   ***                                     \n\n"
                         "   Aligns: %11%\n"
                         "   Coverage (Synthetic): %12%\n"
                         "   Coverage (Genome):    %13%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % "-"
                                            % stats.n_syn
                                            % stats.n_gen
                                            % o.rAnnot
                                            % "??"
                                            % meth2Str()
                                            % "??"
                                            % stats.sBefore
                                            % stats.gBefore
                                            % "??"
                                            % stats.sAfter
                                            % stats.gAfter).str());
    o.writer->close();
}

void RSample::report(const FileName &file, const Options &o)
{
    const auto stats = RSample::stats(file, o);
    
    /*
     * Generating RnaSubsample_summary.stats
     */
    
    generateSummary("RnaSubsample_summary.stats", stats, o);
}