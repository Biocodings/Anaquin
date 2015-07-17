#include "var/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

VAlign::Stats VAlign::analyze(const std::string &file, const Options &options)
{
    VAlign::Stats stats;
    static const auto &s = Standard::instance();

    options.info("Parsing alignment file");

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if (!align.i && (p.i % 1000000) == 0)
        {
            options.wait(std::to_string(p.i));
        }
        
        if (!align.mapped)
        {
            return;
        }
        else if (align.id != s.id)
        {
            stats.n_genome++;
            return;
        }
        
        stats.n_chrT++;

        const Variation *matched;

        if (classify(stats.p.m, align, [&](const Alignment &)
        {
            matched = find(s.v_vars, align, MatchRule::Contains);
            return options.filters.count(matched->id) ? Ignore : matched ? Positive : Negative;
        }))
        {
            stats.c.at(matched->id)++;
        }
    });

    sums(stats.c, stats.p.m.nr);

    /*
     * Generate an abundance plot for the accuracy of quantification, the measured DNA
     * standard abundance (in FPKM) relative to the known concentration (in attamoles/ul)
     * of each DNA standard.
     */

    // The total number of reads aligned
    const auto n = stats.n_chrT;

    // Known concentration for the given mixture
    const auto &m = s.v_seq(options.mix);

    for (const auto &i : stats.c)
    {
        if (!i.second)
        {
            continue;
        }
        
        const auto s = m.at(i.first);
        
        // Compare the FPKM with the known concentration
        const auto known = s.abund();
        
        // Calculate FPKM for the sequin
        const double measured = (std::pow(10, 9) * static_cast<double>(i.second)) / (n * s.l.length());

        stats.z.push_back(i.first);
        stats.x.push_back(log2(known));
        stats.y.push_back(log2(measured));
    }
    
    // Perform a linear regreession
    stats.linear();

    options.info("Calculating LOS");
    
    // Calculate for the sensitivity
    stats.p.s = Expression::analyze(stats.c, s.v_seq(options.mix));

    AnalyzeReporter::stats("var_align.stats", stats.p, stats.c, options.writer);
    
	return stats;
}