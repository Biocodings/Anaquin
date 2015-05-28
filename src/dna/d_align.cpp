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
                stats.c[matched->id]++;
            }
        }
        else
        {
            stats.n_samps++;
        }
    });

    sums(stats.c, stats.p.m.nr);
    //stats.se = Expression::analyze(stats.ce, seqs);

    AnalyzeReporter::report("dna_align.stats", stats, stats.p.m, stats.p.s, stats.c, options.writer);

    std::cout << "Sensitivity: " << stats.p.m.sn() << std::endl;
    std::cout << "Specificity: " << stats.p.m.sp() << std::endl;

	return stats;
}