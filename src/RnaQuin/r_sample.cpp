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
    const auto r = Sampler::subsample(file, stats.p, o);

    stats.after  = r.after;
    stats.before = r.before;
    
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
    
    const auto summary = "-------RnaSubsample Summary Statistics\n\n"
                         "       User generated alignment: %1%\n\n"
                         "-------User alignments (before subsampling)\n\n"
                         "       Synthetic: %2% reads\n"
                         "       Genome:    %3% reads\n\n"
                         "       Method:        %4%\n"
                         "       Normalization: %5%\n\n"
                         "-------User alignments (after subsampling)\n\n"
                         "       Synthetic: %6% reads\n"
                         "       Genome:    %7% reads\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % stats.before.syn
                                            % stats.before.gen
                                            % stats.p
                                            % meth2Str()
                                            % stats.after.syn
                                            % stats.after.gen).str());
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