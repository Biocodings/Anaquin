#include "var/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

VAlign::Stats VAlign::analyze(const std::string &file, const Options &options)
{
    VAlign::Stats stats;
    static const auto &s = Standard::instance();

    options.info("Parsing alignment file");

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
            stats.n_genome++;
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

        stats.x.push_back(log2(known));
        stats.y.push_back(log2(measured));
        stats.z.push_back(i.first);
    }
    
    // Perform a linear regreession
    stats.linear();

    options.info("Calculating LOS");
    
    // Calculate for the sensitivity
    stats.p.s = Expression::analyze(stats.c, s.d_seq(options.mix));

    AnalyzeReporter::stats("var_align.stats", stats.p, stats.c, options.writer);
    
	return stats;
}