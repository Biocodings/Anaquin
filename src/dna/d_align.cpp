#include "dna/d_align.hpp"
#include "stats/expression.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

DAlign::Stats DAlign::analyze(const std::string &file, const Options &options)
{
    DAlign::Stats stats;
    static const auto &s = Standard::instance();

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        const Variation *matched;

        if (align.id == s.id)
        {
            stats.n_chrT++;

            if (classify(stats.p.m, align, [&](const Alignment &)
            {
                matched = findMap(s.d_vars, align, MatchRule::Contains);

                if (options.filters.count(matched->id))
                {
                    return Ignore;
                }
                
                return matched ? Positive : Negative;
            }))
            {
                stats.c.at(matched->id)++;
            }
        }
        else
        {
            stats.n_samps++;
        }
    });

    sums(stats.c, stats.p.m.nr);

    /*
     * Generate an abundance plot for the accuracy of quantification, the measured DNA
     * standard abundance (in RPKM) relative to the known concentration (in attamoles/ul)
     * of each DNA standard.
     */

    // The total number of reads aligned
    const auto n = stats.n_chrT;

    // Known concentration for the given mixture
    const auto &m = s.d_seq(options.mix);

    for (const auto &i : stats.c)
    {
        if (!i.second)
        {
            continue;
        }
        
        const auto s = m.at(i.first);
        
        // Compare the RPKM with the known concentration
        const auto known = s.abund();
        
        // Calculate RPKM for the sequin
        const double measured = (std::pow(10, 9) * static_cast<double>(i.second)) / (n * s.l.length());

        stats.x.push_back(log(known));
        stats.y.push_back(log(measured));
        stats.z.push_back(i.first);
    }
    
    // Perform a linear regreession
    stats.linear();

    // Calculate for the sensitivity
    stats.p.s = 0; // TODO: Fix this!!! Expression::analyze(stats.c, s.d_seq(options.mix));

    AnalyzeReporter::report("dalign.stats", stats.p, stats.c, options.writer);
    
	return stats;
}