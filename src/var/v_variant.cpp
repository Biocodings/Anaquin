#include "var/v_variant.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Spike;

DVariant::Stats DVariant::analyze(const std::string &file, const Options &options)
{
    DVariant::Stats stats;
    const auto &s = Standard::instance();

    ParserVCF::parse(file, [&](const VCFVariant &var, const ParserProgress &)
    {
        Variation match;

        if (classify(stats.m, var, [&](const VCFVariant &)
        {
            // Can we find this variant?
            if (!s.d_vars.count(var.l))
            {
                return Negative;
            }

            match = s.d_vars.at(var.l);

            // Does the variant match with the meta?
            if (match.type != var.type || match.alt != var.alt || match.ref != var.ref)
            {
                return Negative;
            }

            assert(s.d_seqs_bA.count(match.id));
            
            const auto &base = s.d_seqs_bA.at(match.id);
            
            /*
             * Plotting the relative allele frequency that is established by differences
             * in the concentration of reference and variant DNA standards.
             */
            
            // The measured coverage is the number of base calls aligned and used in variant calling
            const auto measured = (double) var.dp_a / (var.dp_r + var.dp_a);

            // The known coverage for allele frequnece
            const auto known = base.alleleFreq();

            stats.x.push_back(known);
            stats.y.push_back(measured);
            stats.z.push_back(match.id);
  
            return Positive;
        }))
        {
            stats.c.at(match.id)++;
        }
    });

    stats.m.nr = s.d_vars.size();
    
    /*
     * Calculate the proportion of genetic variation with alignment coverage
     */
    
    stats.covered = std::accumulate(stats.c.begin(), stats.c.end(), 0,
            [&](int sum, const std::pair<SequinID, Counts> &p)
            {
                return sum + (p.second ? 1 : 0);
            });

    // The proportion of genetic variation with alignment coverage
    stats.covered = stats.covered / s.d_vars.size();

    assert(stats.covered >= 0 && stats.covered <= 1.0);

    // Measure of variant detection independent to sequencing depth or coverage
    stats.efficiency = stats.m.sn() / stats.covered;
    
    /*
     * Created a script for allele frequency
     */

    // Fit a linear model on allele frequency
    stats.linear();
    
    // Create a script for allele frequency
    AnalyzeReporter::script("dna_allele.R", stats.x, stats.y, stats.z, "Allele Frequence", 10.0, options.writer);

    /*
     * Write out results
     */

    const std::string format = "%1%\t%2%\t%3%";

    options.writer->open("dna_variant.stats");
    options.writer->write((boost::format(format) % "sn" % "sp" % "detect").str());
    options.writer->write((boost::format(format) % stats.m.sn()
                                                 % stats.m.sp()
                                                 % stats.covered).str());
    options.writer->write("\n");
    
    for (const auto &p : stats.c)
    {
        options.writer->write((boost::format("%1%\t%2%") % p.first % p.second).str());
    }
    
    options.writer->close();

    return stats;
}