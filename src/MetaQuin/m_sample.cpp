#include "tools/errors.hpp"
#include "MetaQuin/MetaQuin.hpp"
#include "MetaQuin/m_sample.hpp"
#include "writers/sam_writer.hpp"

using namespace Anaquin;

MSample::Stats MSample::analyze(const FileName &file, const Options &o)
{
    A_CHECK(!isnan(o.p), "Sampling probability must not be NAN");
    A_CHECK(o.p > 0 && o.p < 1.0, "Sampling probability must be (0:1)");

    MSample::Stats stats;
    
    o.info(file);
    o.info("Sampling proportion: " + std::to_string(o.p));

    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logInfo(std::to_string(info.p.i));
        }
        
        // Don't count for multiple alignments
        if (x.isPrimary)
        {
            if (isMetaQuin(x.cID))
            {
                stats.before.syn++;
            }
            else
            {
                stats.before.gen++;
            }
        }
    });

    o.info("Alignments mapped to the synthetic community (before subsampling): " + std::to_string(stats.before.syn));
    o.info("Alignments mapped to the genome (before subsampling): " + std::to_string(stats.before.gen));
    
    if (stats.before.syn == 0) { throw std::runtime_error("No alignment found on the in-silico community"); }
    if (stats.before.gen == 0) { throw std::runtime_error("No alignment found on the microbes");   }

    o.info("Calculating the normalization factor");
    
    const auto nTotal = stats.before.gen / (1.0 - o.p);
    A_CHECK(nTotal >= stats.before.gen, "New total is less than number of genomic reads");
    
    // Number of synthetic reads after sampling (eg: 0.10101)
    const auto nSyn = nTotal - stats.before.gen;

    /*
     * Make sure we only derive normalization factor if there're enough synthetic reads.
     */
    
    stats.norm = nSyn < stats.before.syn ? static_cast<Proportion>(nSyn) / stats.before.syn : 1.0;

    o.logInfo("New Total: "     + std::to_string(nTotal));
    o.logInfo("New Synthetic: " + std::to_string(nSyn));
    o.info("Normalization: "    + std::to_string(stats.norm));

    // Perform subsampling
    const auto r = Sampler::sample(file, stats.norm, o, [&](const ChrID &id) { return isMetaQuin(id); });

    stats.after = r.after;

    return stats;
}

static void generateSummary(const FileName &file, const FileName &src, const MSample::Stats &stats, const MSample::Options &o)
{
    o.generate(file);
    
    const auto summary = "-------MetaSubsample Summary Statistics\n\n"
                         "       User generated alignment: %1%\n\n"
                         "-------User alignments (before subsampling)\n\n"
                         "       Synthetic: %2% reads\n"
                         "       Genome:    %3% reads\n"
                         "       Dilution:  %4%\n\n"
                         "       * Dilution specified by the user:\n"
                         "       Fraction: %5%\n\n"
                         "       * Normalization applied in subsampling:\n"
                         "       Normalization: %6%\n\n"
                         "-------User alignments (after subsampling)\n\n"
                         "       Synthetic: %7% reads\n"
                         "       Genome:    %8% reads\n"
                         "       Dilution:  %9%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % src
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

void MSample::report(const FileName &file, const Options &o)
{
    const auto stats = MSample::analyze(file, o);
    
    /*
     * Generating MetaSubsample_summary.stats
     */
    
    generateSummary("MetaSubsample_summary.stats", file, stats, o);
}
