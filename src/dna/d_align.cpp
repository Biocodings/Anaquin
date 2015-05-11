#include <assert.h>
#include "d_align.hpp"
#include "expression.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

DAlignStats DAlign::analyze(const std::string &file, const Options &options)
{
    DAlignStats stats;
    const auto &s = Standard::instance();

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        const BedFeature *matched;
        
        /*
         * Classify at the exon-level
         */
        
        if (classify(stats.me, align, [&](const Alignment &)
                     {
                         matched = find(s.d_exons, align, MatchRule::Contains);
                        
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
            stats.ce[matched->id]++;
        }
    });

    count_ref(stats.ce, stats.me.nr);

    /*
     * Calculate for the LOS
     */
    
    //stats.se = Expression::analyze(stats.ce, seqs);

    AnalyzeReporter::report("dalign_exons.stats", stats.me, stats.se, stats.ce, options.writer);

	return stats;
}