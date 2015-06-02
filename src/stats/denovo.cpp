#include "stats/denovo.hpp"
#include "stats/histogram.h"
#include "parsers/parser_fa.hpp"

using namespace Spike;

DNStats DNAsssembly::stats(const std::string &file)
{
    DNStats stats;
    Histogram h;

    ParserFA::parse(file, [&](const FALine &l, const ParserProgress &)
    {
        Contig c;

        c.id = l.id;
        
        // Sequence of the config
        c.seq = l.seq;
        
        // Length of the contig
        h.insert(c.l = l.seq.length());

        stats.contigs.push_back(c);
    });

    /*
     * This is copied from printContiguityStats() in Histogram.h of the Abyss source code.
     */

    h = h.trimLow(500);
    
    /*
     * Reference: https://github.com/bcgsc/abyss/blob/e58e5a6666e0de0e6bdc15c81fe488f5d83085d1/Common/Histogram.h
     */

    stats.sum  = h.sum();
    stats.N50  = h.n50();
    stats.min  = h.minimum();
    stats.max  = h.maximum();
    stats.mean = h.expectedValue();
    stats.N80  = h.weightedPercentile(1 - 0.8);
    stats.N20  = h.weightedPercentile(1 - 0.2);

    return stats;
}