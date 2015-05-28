#include "expression.hpp"
#include "dna/d_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

DAlignStats DAlign::analyze(const std::string &file, const Options &options)
{
    DAlignStats stats;
    static const auto &s = Standard::instance();

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        const BedFeature *matched;

        if (align.id == s.id)
        {
            stats.n_seqs++;

            if (classify(stats.p.m, align, [&](const Alignment &)
                         {
                             matched = find(s.d_annot, align, MatchRule::Contains);
                             
                             if (!matched)
                             {
                                 return Negative;
                             }
                             else if (options.filters.count(matched->id))
                             {
                                 return Ignore;
                             }
                             
                             return Positive;
                         }))
            {
                stats.c.at(matched->name)++;
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
    const auto n = stats.n_seqs;

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
    stats.p.s = Expression::analyze(stats.c, s.d_seq(options.mix));

    AnalyzeReporter::report("dalign.stats", stats, stats.p.m, stats.p.s, stats.c, options.writer);

    // Generate a R script for plotting relative abundance
    AnalyzeReporter::script("dalign.R", stats, options.writer);

    std::cout << "Sensitivity: " << stats.p.m.sn() << std::endl;
    std::cout << "Specificity: " << stats.p.m.sp() << std::endl;

	return stats;
}