#include <math.h>
#include <iostream>
#include <assert.h>
#include "Biology.hpp"
#include "ParserBED.hpp"
#include "ParserSAM.hpp"
#include "Statistics.hpp"
#include "AlignerAnalyst.hpp"
#include "StandardFactory.hpp"

static bool matchGeneBoundary(const Standard &r, const Alignment &align)
{
    for (auto f: r.fs) // TODO: Double check no intron is here.... // TODO: Check spliced reads
    {
        if (r.l.contains(align.l))
        {
            return true;
        }
    }
    
    return false;
}

static bool matchChromoBoundary(const Standard &r, const Alignment &align)
{
    return r.l.contains(align.l);
}

AlignerStats AlignerAnalyst::spliced(const std::string &file, Sequins s, Reads n)
{
    AlignerStats stats;
    const auto r = StandardFactory::reference();

    ParserBED::parse(file, [&](const BedFeature &f)
    {
        /*
         * The parser simply reports whatever in the file. It doesn't know anything about introns.
         */
        
        Feature t = f;
        
        extractIntrons(f.blocks, [&](BasePair start, BasePair end)
        {
            t.l.set(start, end);
            classify(r.introns, r, t, stats.m);
        });
    });

    return stats;
}

AlignerStats AlignerAnalyst::base(const std::string &file, Sequins s, Reads n)
{
    AlignerStats stats;
    const auto r = StandardFactory::reference();

    /*
     * Calculate the sensitivity and specificity for the experiment. They are statistical measures of the performance
     * of a binary test. Sensitivity (also known true-positive rate) measures the proportion of actual positives which
     * are correctly identified as such (eg: the percentage of reads aligned to the reference correctly). Specificity
     * (also known as true-negative rate) measures the proportion of negatives which are correctly identified as such
     * (eg: the percentage of reads fails to align to the reference incorrectly).
     */
    
    Reads i = 0;
    
    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if (++i > n)
        {
            return false;
        }
        
        // If it's aligned to the reference
        else if (align.id == r.id)
        {
            stats.nr++;
        }
        else
        {
            stats.nq++;
        }
        
		if (align.mapped)
		{
			if (align.id == r.id)
			{
				if (matchChromoBoundary(r, align))
				{
                    // TODO: Rename to feature
					if (matchGeneBoundary(r, align))
					{
						stats.m.tp++;
					}
					else
					{
						stats.m.fp++;
					}
				}
				else
				{
					stats.m.fp++;
				}
			}
			else
			{
				if (matchChromoBoundary(r, align))
				{
					stats.m.fn++;
				}
				else
				{
					stats.m.tn++;
				}
			}
		}

        stats.n++;        
        return true;
    });
    
    assert(stats.n);
    
	stats.pr = static_cast<Percentage>(stats.nr / stats.n);
	stats.pq = static_cast<Percentage>(stats.nq / stats.n);
    stats.dilution = stats.nq ? static_cast<Percentage>(stats.nr / stats.nq) : 1;

	return stats;
}