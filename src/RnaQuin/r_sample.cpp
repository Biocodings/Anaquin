#include "tools/sample.hpp"
#include "RnaQuin/r_sample.hpp"
#include "writers/writer_sam.hpp"

using namespace Anaquin;

RSample::Stats RSample::stats(const FileName &file, const Options &o)
{
    RSample::Stats stats;
    
    o.info(file);

    assert(o.p > 0 && o.p < 1.0);
    o.info("Spike-in proportion: " + std::to_string(o.p));

    /*
     * Calculating sequencing depth for both genomic and synthetic before subsampling
     */

    o.info("[1/3]: Calculating the coverage before subsampling");
    
    ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logInfo(std::to_string(info.p.i));
        }
        
        if (Standard::isSynthetic(x.cID))
        {
            stats.before.syn++;
        }
        else
        {
            stats.before.gen++;
        }
    });

    o.info("Alignments mapped to the in silico (before subsampling): " + std::to_string(stats.before.syn));
    o.info("Alignments mapped to the genome (before subsampling): " + std::to_string(stats.before.gen));
    
    if (stats.before.syn == 0) { throw std::runtime_error("No alignment found on the in silico chromosome"); }
    if (stats.before.gen == 0) { throw std::runtime_error("No alignment found on the genome");   }

    /*
     * Calculating the subsamping fraction. Eg: if we have 10m reads to the genome and 5m reads to the
     * in-silico chromsome. The specified fraction is 1%.
     *
     *   New total is: 10/0.99 = 10.10101 => Synthetic will have: 0.10101 reads.
     *   We should sample for 0.10101/5 = 0.020202.
     */
    
    o.info("[2/3]: Calculating the normalization factor");
    
    const auto nTotal = stats.before.gen / (1.0 - o.p);
    o.logInfo("New Total: " + std::to_string(nTotal));
    
    assert(nTotal > stats.before.gen);
    
    const auto nSyn = nTotal - stats.before.gen;
    o.logInfo("New Synthetic: " + std::to_string(nSyn));

    stats.norm = nSyn < stats.before.syn ? static_cast<Proportion>(nSyn) / stats.before.syn : 1.0;
    o.info("Normalization: " + std::to_string(stats.norm));
    
    // Perform subsampling
    const auto r = Sampler::subsample(file, stats.norm, o);

    stats.after = r.after;

    return stats;
}

static void generateSummary(const FileName &file, const RSample::Stats &stats, const RSample::Options &o)
{
    o.generate(file);
    
    const auto summary = "-------RnaSubsample Summary Statistics\n\n"
                         "       User generated alignment: %1%\n\n"
                         "-------User alignments (before subsampling)\n\n"
                         "       Synthetic: %2% reads\n"
                         "       Genome:    %3% reads\n"
                         "       Dilution:  %4%\n\n"
                         "       * Fraction of dilution specified:\n"
                         "       Fraction: %5%\n\n"
                         "       * Normalization applied in subsampling:\n"
                         "       Normalization: %6%\n\n"
                         "-------User alignments (after subsampling)\n\n"
                         "       Synthetic: %7% reads\n"
                         "       Genome:    %8% reads\n"
                         "       Dilution:  %9%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % stats.before.syn
                                            % stats.before.gen
                                            % stats.before.dilut()
                                            % std::to_string(o.p)
                                            % stats.norm
                                            % stats.after.syn
                                            % stats.after.gen
                                            % stats.after.dilut()).str());
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