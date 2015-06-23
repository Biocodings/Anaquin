#include "dna/d_variant.hpp"
#include "stats/expression.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Spike;

DVariant::Stats DVariant::analyze(const std::string &file, const Options &options)
{
    DVariant::Stats stats;
    const auto &s = Standard::instance();

    ParserVCF::parse(file, [&](const VCFVariant &var, const ParserProgress &)
    {
        Variation known;

        if (classify(stats.m, var, [&](const VCFVariant &)
        {
            // Can we find this variant?
            if (!s.d_vars.count(var.l))
            {
                return Negative;
            }

            known = s.d_vars.at(var.l);

            // Does the variant match with the meta?
            if (known.type != var.type || known.alt != var.alt || known.ref != var.ref)
            {
                return Negative;
            }

            return Positive;
        }))
        {
            stats.c.at(known.id)++;
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