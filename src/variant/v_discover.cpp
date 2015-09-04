#include "variant/v_discover.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

VDiscover::Stats VDiscover::analyze(const std::string &file, const Options &options)
{
    VDiscover::Stats stats;
    const auto &s = Standard::instance();

    options.info("Parsing VCF file");

    ParserVCF::parse(file, [&](const VCFVariant &var, const ParserProgress &)
    {
        Variation match;

        if (classify(stats.m, var, [&](const VCFVariant &)
        {
            // Can we find this variant?
            if (!s.v_vars.count(var.l))
            {
                return Negative;
            }

            match = s.v_vars.at(var.l);

            // Does the variant match with the meta?
            if (match.type != var.type || match.alt != var.alt || match.ref != var.ref)
            {
                return Negative;
            }

            assert(s.bases_1.count(match.id));
            
            return Positive;
        }))
        {
            stats.h.at(match)++;
        }
    });

    options.info("Generating statistics");

    stats.m.nr = s.v_vars.size();

    /*
     * Calculate the proportion of genetic variation with alignment coverage
     */
    
    stats.covered = std::accumulate(stats.h.begin(), stats.h.end(), 0,
            [&](int sum, const std::pair<Locus, Counts> &p)
            {
                return sum + (p.second ? 1 : 0);
            });

    // The proportion of genetic variation with alignment coverage
    stats.covered = stats.covered / s.v_vars.size();

    assert(stats.covered >= 0 && stats.covered <= 1.0);

    // Measure of variant detection independent to sequencing depth or coverage
    stats.efficiency = stats.m.sn() / stats.covered;
    
    /*
     * Generate summary statistics
     */

    const std::string format = "%1%\t%2%\t%3%";
    
    options.writer->open("VarDiscover_summary.stats");
    options.writer->write((boost::format(format) % "sn" % "sp" % "detect").str());
    options.writer->write((boost::format(format) % stats.m.sn()
                           % stats.m.sp()
                           % stats.covered).str());
    
    for (const auto &p : stats.h)
    {
        options.writer->write((boost::format("%1%-%2%\t%3%") % p.first.id % p.first.l.start % p.second).str());
    }
    
    options.writer->close();

    return stats;
}